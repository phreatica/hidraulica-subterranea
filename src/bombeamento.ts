/**
 * Ensaios de bombeamento (vazão constante) em aquífero confinado. Unidades SI: m, s, m³/s; T em m²/s.
 *
 * Theis (1935):          s = Q / (4πT) · W(u),   u = r²S / (4Tt),   W(u) = E₁(u) (integral exponencial)
 * Cooper-Jacob (1946):   para u pequeno (u < 0,01 é a recomendação usual),
 *                        T = 2,303·Q / (4π·Δs)   (Δs = rebaixamento por ciclo log de t)
 *                        S = 2,25·T·t₀ / r²      (t₀ = tempo em que a reta corta s = 0)
 * Recuperação de Theis:  s' = 2,303·Q/(4πT) · log10(t/t′) → T = 2,303·Q / (4π·Δs′)
 * Hantush-Jacob (1955), aquífero semiconfinado sem armazenamento no aquitarde:
 *                        s = Q/(4πT) · W(u, r/B),  W(u, β) = ∫ᵤ^∞ (1/y)·exp(−y − β²/(4y)) dy
 *
 * Referências: Theis, C. V. (1935) Trans. AGU 16:519–524; Cooper, H. H.; Jacob, C. E. (1946) Trans. AGU 27:526–534;
 * Hantush, M. S.; Jacob, C. E. (1955) Trans. AGU 36:95–100; Kruseman, G. P.; De Ridder, N. A. (1990) Analysis and
 * Evaluation of Pumping Test Data, ILRI Publ. 47.
 */
import { regressaoLinear, type Reta } from "./estatistica";

const EULER = 0.5772156649015329;

/** Integral exponencial E₁(x) = W(u) de Theis. Série para x ≤ 1; fração contínua (Lentz) para x > 1. */
export function E1(x: number): number {
  if (!(x > 0)) throw new RangeError("E1 está definida para x > 0.");
  if (x <= 1) {
    let soma = 0;
    let termo = 1;
    for (let k = 1; k < 200; k++) {
      termo *= -x / k;
      const parcela = -termo / k;
      soma += parcela;
      if (Math.abs(parcela) < 1e-17 * Math.abs(soma)) break;
    }
    return -EULER - Math.log(x) + soma;
  }
  // fração contínua: E1(x) = e^{-x} · 1/(x+1−1²/(x+3−2²/(x+5−…)))
  const TINY = 1e-300;
  let b = x + 1;
  let c = 1 / TINY;
  let d = 1 / b;
  let h = d;
  for (let i = 1; i < 300; i++) {
    const a = -i * i;
    b += 2;
    d = 1 / (a * d + b);
    c = b + a / c;
    const del = c * d;
    h *= del;
    if (Math.abs(del - 1) < 1e-16) break;
  }
  return h * Math.exp(-x);
}

/** Função de poço de Theis W(u). */
export const funcaoPocoTheis = E1;

/** Rebaixamento de Theis s (m) para Q (m³/s), T (m²/s), S, r (m) e t (s). */
export function rebaixamentoTheis(p: { Q_m3_s: number; T_m2_s: number; S: number; r_m: number; t_s: number }): number {
  const { Q_m3_s, T_m2_s, S, r_m, t_s } = p;
  const u = (r_m * r_m * S) / (4 * T_m2_s * t_s);
  return (Q_m3_s / (4 * Math.PI * T_m2_s)) * E1(u);
}

/** u de Theis. */
export function uTheis(p: { T_m2_s: number; S: number; r_m: number; t_s: number }): number {
  return (p.r_m * p.r_m * p.S) / (4 * p.T_m2_s * p.t_s);
}

export interface LeituraBombeamento {
  t_s: number;
  s_m: number; // rebaixamento
}

export interface ResultadoCooperJacob {
  T_m2_s: number;
  S: number;
  reta: Reta; // s = a + b·log10(t)
  deltaS_m: number; // rebaixamento por ciclo log
  t0_s: number;
  uMaximo: number; // maior u entre os pontos usados (validade da aproximação)
  avisos: string[];
}

/** Cooper-Jacob: reta s × log10(t) nos pontos com t ≥ tMinimo. */
export function cooperJacob(leituras: LeituraBombeamento[], p: { Q_m3_s: number; r_m: number; tMinimo_s?: number }): ResultadoCooperJacob {
  const pts = leituras.filter((l) => l.t_s > 0 && Number.isFinite(l.s_m) && l.t_s >= (p.tMinimo_s ?? 0));
  if (pts.length < 2) throw new RangeError("São necessárias pelo menos 2 leituras no trecho escolhido.");
  const reta = regressaoLinear(pts.map((l) => Math.log10(l.t_s)), pts.map((l) => l.s_m));
  const deltaS = reta.b;
  if (!(deltaS > 0)) throw new RangeError("O rebaixamento não cresce com o tempo no trecho escolhido.");
  const T = (2.303 * p.Q_m3_s) / (4 * Math.PI * deltaS);
  const t0 = Math.pow(10, -reta.a / reta.b);
  const S = (2.25 * T * t0) / (p.r_m * p.r_m);
  const tMin = Math.min(...pts.map((l) => l.t_s));
  const uMax = uTheis({ T_m2_s: T, S, r_m: p.r_m, t_s: tMin });
  const avisos: string[] = [];
  if (uMax > 0.01) avisos.push(`u = ${uMax.toPrecision(2)} no primeiro ponto usado: acima de 0,01, a aproximação de Cooper-Jacob perde precisão. Use só os tempos maiores.`);
  return { T_m2_s: T, S, reta, deltaS_m: deltaS, t0_s: t0, uMaximo: uMax, avisos };
}

/** Recuperação de Theis: reta s′ × log10(t/t′), onde t = tempo desde o início do bombeamento e t′ = desde a parada. */
export function recuperacaoTheis(leituras: { t_s: number; tLinha_s: number; sResidual_m: number }[], Q_m3_s: number): { T_m2_s: number; reta: Reta; deltaS_m: number } {
  const pts = leituras.filter((l) => l.tLinha_s > 0 && l.t_s > l.tLinha_s);
  if (pts.length < 2) throw new RangeError("São necessárias pelo menos 2 leituras de recuperação.");
  const reta = regressaoLinear(pts.map((l) => Math.log10(l.t_s / l.tLinha_s)), pts.map((l) => l.sResidual_m));
  if (!(reta.b > 0)) throw new RangeError("O rebaixamento residual não cresce com t/t′: verifique os dados.");
  return { T_m2_s: (2.303 * Q_m3_s) / (4 * Math.PI * reta.b), reta, deltaS_m: reta.b };
}

/** Ajuste de Theis por mínimos quadrados em log(T) e log(S) (Gauss-Newton com derivadas numéricas). */
export function ajustarTheis(
  leituras: LeituraBombeamento[],
  p: { Q_m3_s: number; r_m: number; T0_m2_s?: number; S0?: number },
): { T_m2_s: number; S: number; rmse_m: number; iteracoes: number } {
  const pts = leituras.filter((l) => l.t_s > 0 && l.s_m > 0);
  if (pts.length < 3) throw new RangeError("São necessárias pelo menos 3 leituras com rebaixamento positivo.");
  let x = [Math.log(p.T0_m2_s ?? 1e-3), Math.log(p.S0 ?? 1e-4)];
  if (p.T0_m2_s === undefined) {
    try {
      const cj = cooperJacob(pts, { Q_m3_s: p.Q_m3_s, r_m: p.r_m, tMinimo_s: pts[Math.floor(pts.length / 2)].t_s });
      if (cj.T_m2_s > 0 && cj.S > 0) x = [Math.log(cj.T_m2_s), Math.log(cj.S)];
    } catch {
      /* mantém o chute inicial */
    }
  }
  const residuos = (v: number[]) =>
    pts.map((l) => rebaixamentoTheis({ Q_m3_s: p.Q_m3_s, T_m2_s: Math.exp(v[0]), S: Math.exp(v[1]), r_m: p.r_m, t_s: l.t_s }) - l.s_m);
  const custo = (v: number[]) => residuos(v).reduce((s, r) => s + r * r, 0);
  let it = 0;
  for (; it < 100; it++) {
    const r0 = residuos(x);
    const h = 1e-6;
    const J = pts.map(() => [0, 0]);
    for (let k = 0; k < 2; k++) {
      const xp = x.slice();
      xp[k] += h;
      const rp = residuos(xp);
      for (let i = 0; i < pts.length; i++) J[i][k] = (rp[i] - r0[i]) / h;
    }
    let a11 = 0, a12 = 0, a22 = 0, g1 = 0, g2 = 0;
    for (let i = 0; i < pts.length; i++) {
      a11 += J[i][0] * J[i][0]; a12 += J[i][0] * J[i][1]; a22 += J[i][1] * J[i][1];
      g1 += J[i][0] * r0[i]; g2 += J[i][1] * r0[i];
    }
    const det = a11 * a22 - a12 * a12;
    if (Math.abs(det) < 1e-30) break;
    const d0 = -(a22 * g1 - a12 * g2) / det;
    const d1 = -(-a12 * g1 + a11 * g2) / det;
    // busca linear simples para garantir descida
    let passo = 1;
    const c0 = custo(x);
    let novo = [x[0] + d0, x[1] + d1];
    while (custo(novo) > c0 && passo > 1e-4) {
      passo /= 2;
      novo = [x[0] + passo * d0, x[1] + passo * d1];
    }
    x = novo;
    if (Math.abs(passo * d0) < 1e-10 && Math.abs(passo * d1) < 1e-10) break;
  }
  const rmse = Math.sqrt(custo(x) / pts.length);
  return { T_m2_s: Math.exp(x[0]), S: Math.exp(x[1]), rmse_m: rmse, iteracoes: it };
}

/**
 * Função de poço de Hantush-Jacob W(u, β), β = r/B, por integração numérica (Simpson em ln y).
 * Para β = 0 coincide com W(u) de Theis.
 */
export function funcaoPocoHantush(u: number, beta: number): number {
  if (!(u > 0)) throw new RangeError("u deve ser maior que zero.");
  if (!(beta >= 0)) throw new RangeError("β não pode ser negativo.");
  if (beta === 0) return E1(u);
  // substituição y = e^z: ∫ exp(−y − β²/(4y)) dz, de z = ln u até ln(u) + ~ (y grande)
  const zIni = Math.log(u);
  const zFim = Math.log(Math.max(u, 1) * 60 + 60);
  const n = 4000; // par
  const hz = (zFim - zIni) / n;
  const f = (z: number) => {
    const y = Math.exp(z);
    return Math.exp(-y - (beta * beta) / (4 * y));
  };
  let soma = f(zIni) + f(zFim);
  for (let i = 1; i < n; i++) soma += (i % 2 ? 4 : 2) * f(zIni + i * hz);
  return (soma * hz) / 3;
}

/** Rebaixamento de Hantush-Jacob (aquífero semiconfinado). B = fator de drenança (m). */
export function rebaixamentoHantush(p: { Q_m3_s: number; T_m2_s: number; S: number; r_m: number; t_s: number; B_m: number }): number {
  const u = uTheis(p);
  return (p.Q_m3_s / (4 * Math.PI * p.T_m2_s)) * funcaoPocoHantush(u, p.r_m / p.B_m);
}
