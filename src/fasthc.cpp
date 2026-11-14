// Gaussian-BIC hill climbing from sufficient statistics.
//
// Everything is computed from the p x p correlation matrix R. For node j with
// parent set P (|P| = k) let A = R[P, P] = L L' (Cholesky) and
// B = L^{-1} R[P, ] (k x p). Then, in correlation units,
//
//   rss_j(P)          = 1 - ||B[, j]||^2
//   partial var of i  : d_i = 1 - ||B[, i]||^2
//   partial cov (i,j) : c_i = R[i, j] - B[, i]' B[, j]
//   add i  to P       : rss_j(P + i) = rss_j(P) - c_i^2 / d_i
//   drop m from P     : rss_j(P - m) = rss_j(P) + beta_m^2 / (A^{-1})_mm,
//                       beta = A^{-1} R[P, j],  diag(A^{-1}) = colSums(L^{-1}^2)
//
// so a whole column of the score-delta cache costs O(k^3 + p k^2), with no
// dependence on n. Moves, ordering and tolerances mirror bnlearn::hc().

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <climits>

using namespace Rcpp;

namespace {

const double FHC_TOL = 1.490116119384765625e-08;   // sqrt(DBL_EPSILON), as bnlearn


struct Problem {
  const arma::mat& R;
  int p;
  double n, logn;
  bool mle;
  double css;
  std::vector<double> cssv;
  Problem(const arma::mat& R_, double n_, bool mle_)
    : R(R_), p(R_.n_rows), n(n_), logn(std::log(n_)), mle(mle_), css(1) {}

  // node score up to node-specific constants (they cancel in all deltas).
  double score(double rss, int k) const {
    double dfs = mle ? n : n - k - 1.0;
    if (!(dfs > 0) || !(rss * css / dfs >= DBL_EPSILON)) return R_NegInf;   // bnlearn: sd < sqrt(eps)
    if (mle)
      return -0.5 * n * std::log(rss) - 0.5 * logn * k;
    double df = n - k - 1.0;
    if (df <= 0) return R_NegInf;
    return -0.5 * n * std::log(rss) + 0.5 * n * std::log(df) - 0.5 * df -
           0.5 * logn * k;
  }
};

// recompute rss of node j and its column of the cache.
double fill_column(Problem pr, int j, const std::vector<int>& pa,
                   const std::vector<char>& bl, arma::mat& cache) {

  const int p = pr.p, k = pa.size();
  pr.css = pr.cssv[j];
  const arma::mat& R = pr.R;
  double rss;

  if (k == 0) {

    rss = 1.0;
    double s0 = pr.score(rss, 0);
    for (int i = 0; i < p; i++) {
      if (i == j || bl[i + j * p]) continue;
      double rij = R(i, j);
      cache(i, j) = pr.score(1.0 - rij * rij, 1) - s0;
    }
    return rss;

  }

  arma::uvec P(k);
  for (int m = 0; m < k; m++) P(m) = pa[m];

  arma::mat A = R.submat(P, P);
  arma::mat L;
  if (!arma::chol(L, A, "lower")) {
    // parents are (numerically) collinear; nothing sensible can be added.
    for (int i = 0; i < p; i++) if (i != j) cache(i, j) = R_NegInf;
    return 0.0;
  }

  arma::mat Linv = arma::solve(arma::trimatl(L), arma::eye(k, k));
  arma::mat B = Linv * R.rows(P);                  // k x p
  arma::rowvec nb = arma::sum(arma::square(B), 0); // ||B[, i]||^2
  arma::vec bj = B.col(j);
  rss = 1.0 - nb(j);

  arma::rowvec cpart = R.row(j) - bj.t() * B;      // partial covariances
  arma::vec beta = Linv.t() * bj;                  // A^{-1} R[P, j]
  arma::rowvec ainv = arma::sum(arma::square(Linv), 0);

  double s0 = pr.score(rss, k);

  std::vector<int> pos(p, -1);
  for (int m = 0; m < k; m++) pos[pa[m]] = m;

  for (int i = 0; i < p; i++) {

    if (i == j || bl[i + j * p]) continue;

    if (pos[i] >= 0) {
      int m = pos[i];
      double rnew = rss + beta(m) * beta(m) / ainv(m);
      cache(i, j) = pr.score(rnew, k - 1) - s0;
    }
    else {
      double d = 1.0 - nb(i);
      if (d < 1e-15) { cache(i, j) = R_NegInf; continue; }
      double rnew = rss - cpart(i) * cpart(i) / d;
      cache(i, j) = pr.score(rnew, k + 1) - s0;
    }

  }

  return rss;

}

// is there a directed path from -> to? if notdirect, ignore the arc from -> to.
bool has_path(int from, int to, const std::vector<std::vector<int> >& ch,
              std::vector<int>& stack, std::vector<char>& seen, bool notdirect) {

  std::fill(seen.begin(), seen.end(), 0);
  stack.clear();
  stack.push_back(from);
  seen[from] = 1;

  while (!stack.empty()) {
    int v = stack.back(); stack.pop_back();
    for (size_t c = 0; c < ch[v].size(); c++) {
      int w = ch[v][c];
      if (notdirect && v == from && w == to) continue;
      if (w == to) return true;
      if (!seen[w]) { seen[w] = 1; stack.push_back(w); }
    }
  }

  return false;

}

void erase_val(std::vector<int>& v, int x) {
  for (size_t m = 0; m < v.size(); m++)
    if (v[m] == x) { v.erase(v.begin() + m); return; }
}

} // namespace

// [[Rcpp::export(name = ".fasthc_cpp")]]
List fasthc_cpp(const arma::mat& R, NumericVector css, double n, IntegerMatrix blmat,
                IntegerMatrix wlmat, double maxp, double maxiter, bool mle) {

  const int p = R.n_rows;
  Problem pr(R, n, mle);
  pr.cssv.assign(css.begin(), css.end());

  int mp = (std::isfinite(maxp) && maxp < p) ? (int)maxp : p;
  double mi = std::isfinite(maxiter) ? maxiter : R_PosInf;

  std::vector<char> am(p * p, 0), bl(p * p, 0), wl(p * p, 0);
  for (int c = 0; c < p * p; c++) { bl[c] = blmat[c] != 0; wl[c] = wlmat[c] != 0; }

  std::vector<std::vector<int> > pa(p), ch(p);
  std::vector<int> stack; stack.reserve(p);
  std::vector<char> seen(p);

  // start from the whitelisted arcs.
  for (int i = 0; i < p; i++)
    for (int j = 0; j < p; j++)
      if (wl[i + j * p]) {
        am[i + j * p] = 1; pa[j].push_back(i); ch[i].push_back(j);
      }
  for (int i = 0; i < p; i++)
    for (size_t c = 0; c < ch[i].size(); c++)
      if (has_path(ch[i][c], i, ch, stack, seen, false))
        stop("the whitelist contains a directed cycle.");

  arma::mat cache(p, p, arma::fill::zeros);
  std::vector<double> rss(p);
  for (int j = 0; j < p; j++) rss[j] = fill_column(pr, j, pa[j], bl, cache);

  double iter = 0;
  int nadd = 0, ndrop = 0, nrev = 0;

  while (iter < mi) {

    double best = 0, t;
    int op = 0, from = -1, to = -1;

    // 1) arc additions.
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        if (i == j || am[i + j * p] || am[j + i * p] || bl[i + j * p]) continue;
        if ((int)pa[j].size() >= mp) continue;
        t = cache(i, j);
        if (t - best > FHC_TOL) {
          if (has_path(j, i, ch, stack, seen, false)) continue;
          best = t; op = 1; from = i; to = j;
        }
      }
    }

    // 2) arc deletions.
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        if (!am[i + j * p] || wl[i + j * p]) continue;
        t = cache(i, j);
        if (t - best > FHC_TOL) { best = t; op = 2; from = i; to = j; }
      }
    }

    // 3) arc reversals.
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        if (!am[i + j * p] || bl[j + i * p]) continue;
        if ((int)pa[i].size() >= mp) continue;
        t = cache(i, j) + cache(j, i);
        if (std::fabs(t) < FHC_TOL) t = 0;
        if (t - best > FHC_TOL) {
          if (has_path(i, j, ch, stack, seen, true)) continue;
          best = t; op = 3; from = i; to = j;
        }
      }
    }

    if (op == 0) break;

    if (op == 1) {
      am[from + to * p] = 1; pa[to].push_back(from); ch[from].push_back(to);
      rss[to] = fill_column(pr, to, pa[to], bl, cache);
      nadd++;
    }
    else if (op == 2) {
      am[from + to * p] = 0; erase_val(pa[to], from); erase_val(ch[from], to);
      rss[to] = fill_column(pr, to, pa[to], bl, cache);
      ndrop++;
    }
    else {
      am[from + to * p] = 0; erase_val(pa[to], from); erase_val(ch[from], to);
      am[to + from * p] = 1; pa[from].push_back(to); ch[to].push_back(from);
      rss[from] = fill_column(pr, from, pa[from], bl, cache);
      rss[to] = fill_column(pr, to, pa[to], bl, cache);
      nrev++;
    }

    iter++;

  }

  IntegerMatrix amat(p, p);
  IntegerVector npar(p);
  for (int c = 0; c < p * p; c++) amat[c] = am[c];
  for (int j = 0; j < p; j++) npar[j] = pa[j].size();

  return List::create(_["amat"] = amat, _["rss"] = wrap(rss),
                      _["nparents"] = npar, _["iterations"] = iter,
                      _["moves"] = IntegerVector::create(nadd, ndrop, nrev));

}
