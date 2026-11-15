/**
 * Ensaio de carga ou descarga instantânea (slug test) em poço ou piezômetro.
 * Unidades SI: metros e segundos; K em m/s.
 *
 * Hvorslev (1951), filtro com Le/R > 8:
 *   K = r_c² · ln(Le / R) / (2 · Le · T₀)
 *   T₀ = tempo em que h/h₀ = 0,37 (1/e), obtido da reta ln(h/h₀) × t.
 *
 * Bouwer & Rice (1976), aquífero livre:
 *   K = r_c² · ln(R_e / r_w) / (2 · Le) · (1/t) · ln(y₀ / y_t)
 *   ln(R_e/r_w) = [ 1,1 / ln(Lw/r_w) + (A + B·ln((H − Lw)/r_w)) / (Le/r_w) ]⁻¹          (Lw < H)
 *   ln(R_e/r_w) = [ 1,1 / ln(Lw/r_w) + C / (Le/r_w) ]⁻¹                                 (Lw = H, poço completo)
 *   ln((H − Lw)/r_w) limitado a 6, como no método original.
 *   A, B e C vêm das curvas de Bouwer & Rice (1976); aqui pela aproximação polinomial de
 *   Van Rooy (1988), válida para Le/r_w até cerca de 2 000.
 *
 * Variáveis: r_c = raio do revestimento onde o nível varia; R e r_w = raio do furo/filtro (com pré-filtro);
 * Le = comprimento do filtro; Lw = distância do N.A. estático ao fundo do filtro; H = espessura saturada.
 *
 * Referências: Hvorslev, M. J. (1951) Time lag and soil permeability in ground-water observations, USACE Bull. 36.
 * Bouwer, H.; Rice, R. C. (1976) Water Resources Research 12(3):423–428. Bouwer, H. (1989) Ground Water 27(3):304–309.
 */
import { regressaoLinear, type Reta } from "./estatistica";

export interface LeituraSlug {
  t_s: number; // tempo desde o início (s)
  h_m: number; // deslocamento em relação ao N.A. estático (m), positivo ou negativo
}

export interface AjusteSlug {
  reta: Reta; // ln(h/h₀) = a + b·t
  h0_m: number;
  pontosUsados: number;
}

/** Ajusta ln(h/h₀) × t (com h₀ = primeiro deslocamento), opcionalmente só no trecho de h/h₀ entre `faixa`. */
export function ajustarRecuperacao(leituras: LeituraSlug[], faixa: [number, number] = [0.15, 1]): AjusteSlug {
  const validas = leituras.filter((l) => Number.isFinite(l.t_s) && Number.isFinite(l.h_m) && l.h_m !== 0);
  if (validas.length < 3) throw new RangeError("São necessárias pelo menos 3 leituras com deslocamento diferente de zero.");
  const h0 = validas[0].h_m;
  const pts = validas
    .map((l) => ({ t: l.t_s, r: l.h_m / h0 }))
    .filter((p) => p.r > 0 && p.r >= faixa[0] && p.r <= faixa[1]);
  if (pts.length < 2) throw new RangeError("Poucos pontos na faixa escolhida de h/h₀.");
  const reta = regressaoLinear(pts.map((p) => p.t), pts.map((p) => Math.log(p.r)));
  if (!(reta.b < 0)) throw new RangeError("O deslocamento não diminui com o tempo: verifique os dados.");
  return { reta, h0_m: Math.abs(h0), pontosUsados: pts.length };
}

/** Tempo básico de Hvorslev T₀ (s): tempo em que a reta ajustada passa por h/h₀ = 1/e. */
export function tempoBasicoHvorslev(ajuste: AjusteSlug): number {
  return (-1 - ajuste.reta.a) / ajuste.reta.b;
}

/** K de Hvorslev (m/s). Exige Le/R > 8 (aviso no retorno quando não atende). */
export function hvorslev(p: { rc_m: number; R_m: number; Le_m: number; T0_s: number }): { K_m_s: number; avisos: string[] } {
  const { rc_m, R_m, Le_m, T0_s } = p;
  for (const [n, v] of Object.entries(p)) if (!(v > 0)) throw new RangeError(`${n} deve ser maior que zero.`);
  const avisos: string[] = [];
  if (Le_m / R_m <= 8) avisos.push("Le/R ≤ 8: a fórmula simplificada de Hvorslev supõe Le/R > 8.");
  return { K_m_s: (rc_m * rc_m * Math.log(Le_m / R_m)) / (2 * Le_m * T0_s), avisos };
}

/** Coeficientes A, B e C de Bouwer & Rice pela aproximação polinomial de Van Rooy (1988); x = Le/r_w. */
export function coeficientesBouwerRice(x: number): { A: number; B: number; C: number } {
  if (!(x > 0)) throw new RangeError("Le/r_w deve ser maior que zero.");
  const A = 1.472 + 3.537e-2 * x - 8.148e-5 * x ** 2 + 1.028e-7 * x ** 3 - 6.484e-11 * x ** 4 + 1.573e-14 * x ** 5;
  const B = 0.2372 + 5.151e-3 * x - 2.682e-6 * x ** 2 - 3.491e-10 * x ** 3 + 4.738e-13 * x ** 4;
  const C = 0.792 + 3.993e-2 * x - 5.743e-5 * x ** 2 + 3.858e-8 * x ** 3 - 9.659e-12 * x ** 4;
  return { A, B, C };
}

/** ln(R_e/r_w) de Bouwer & Rice. Lw = H indica poço que penetra todo o aquífero. */
export function lnReRw(p: { rw_m: number; Le_m: number; Lw_m: number; H_m: number }): { valor: number; avisos: string[] } {
  const { rw_m, Le_m, Lw_m, H_m } = p;
  for (const [n, v] of Object.entries(p)) if (!(v > 0)) throw new RangeError(`${n} deve ser maior que zero.`);
  if (Lw_m > H_m + 1e-9) throw new RangeError("Lw não pode ser maior que H.");
  const avisos: string[] = [];
  const x = Le_m / rw_m;
  if (x > 2000) avisos.push("Le/r_w acima de 2 000: fora da faixa da aproximação de Van Rooy.");
  const { A, B, C } = coeficientesBouwerRice(x);
  const termo1 = 1.1 / Math.log(Lw_m / rw_m);
  let termo2: number;
  if (Math.abs(H_m - Lw_m) < 1e-9) {
    termo2 = C / x;
  } else {
    const lnDH = Math.min(Math.log((H_m - Lw_m) / rw_m), 6);
    termo2 = (A + B * lnDH) / x;
  }
  return { valor: 1 / (termo1 + termo2), avisos };
}

/** K de Bouwer & Rice (m/s) a partir da inclinação da reta ln(y) × t (= −(1/t)·ln(y₀/y_t)). */
export function bouwerRice(p: { rc_m: number; rw_m: number; Le_m: number; Lw_m: number; H_m: number; inclinacao_por_s: number }): {
  K_m_s: number;
  lnReRw: number;
  avisos: string[];
} {
  const { rc_m, inclinacao_por_s, ...geo } = p;
  if (!(rc_m > 0)) throw new RangeError("rc_m deve ser maior que zero.");
  if (!(inclinacao_por_s < 0)) throw new RangeError("A inclinação de ln(y) × t deve ser negativa (nível voltando ao estático).");
  const ln = lnReRw(geo);
  const K = ((rc_m * rc_m * ln.valor) / (2 * geo.Le_m)) * -inclinacao_por_s;
  return { K_m_s: K, lnReRw: ln.valor, avisos: ln.avisos };
}
