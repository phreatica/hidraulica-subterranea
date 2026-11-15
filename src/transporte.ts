/**
 * Fluxo e transporte advectivo em meio poroso (unidades SI: m, s, kg).
 *
 * Lei de Darcy:            q = K · i                (fluxo específico, m/s)
 * Velocidade linear média: v = K · i / nₑ           (m/s)
 * Tempo de trânsito:       t = L / v                (s)
 * Fator de retardamento:   R = 1 + ρ_b · K_d / n    (adimensional)
 * Coeficiente de partição: K_d = K_oc · f_oc        (m³/kg)
 * Velocidade do soluto:    v_c = v / R              (m/s)
 *
 * Referências: Fetter, C. W. Applied Hydrogeology; Freeze & Cherry (1979), Groundwater.
 */

function positivo(nome: string, v: number): void {
  if (!Number.isFinite(v) || v <= 0) throw new RangeError(`${nome} deve ser um número maior que zero (recebido: ${v}).`);
}

/** Gradiente hidráulico entre dois pontos: (h₁ − h₂) / L. Positivo quando a carga cai de 1 para 2. */
export function gradienteHidraulico(carga1_m: number, carga2_m: number, distancia_m: number): number {
  positivo("distancia_m", distancia_m);
  return (carga1_m - carga2_m) / distancia_m;
}

/** Fluxo específico de Darcy q = K·i (m/s). */
export function fluxoDarcy(K_m_s: number, gradiente: number): number {
  positivo("K_m_s", K_m_s);
  return K_m_s * gradiente;
}

/** Velocidade linear média (de infiltração) v = K·i/nₑ (m/s). */
export function velocidadeLinear(K_m_s: number, gradiente: number, porosidadeEfetiva: number): number {
  positivo("K_m_s", K_m_s);
  if (!(porosidadeEfetiva > 0 && porosidadeEfetiva <= 1)) throw new RangeError("porosidadeEfetiva deve estar entre 0 e 1.");
  return (K_m_s * gradiente) / porosidadeEfetiva;
}

/** Tempo de trânsito t = L / v (s). Usa o módulo da velocidade. */
export function tempoTransito(distancia_m: number, velocidade_m_s: number): number {
  positivo("distancia_m", distancia_m);
  const v = Math.abs(velocidade_m_s);
  positivo("velocidade", v);
  return distancia_m / v;
}

/** K_d = K_oc · f_oc. K_oc em L/kg (como nas tabelas) → K_d em L/kg. */
export function coeficienteParticao(Koc_L_kg: number, foc: number): number {
  positivo("Koc_L_kg", Koc_L_kg);
  if (!(foc >= 0 && foc <= 1)) throw new RangeError("foc deve estar entre 0 e 1 (fração, não %).");
  return Koc_L_kg * foc;
}

/**
 * Fator de retardamento R = 1 + ρ_b·K_d/n.
 * ρ_b em kg/L (= g/cm³, ex.: 1,6) e K_d em L/kg: o produto é adimensional.
 */
export function fatorRetardamento(densidadeAparente_kg_L: number, Kd_L_kg: number, porosidadeTotal: number): number {
  positivo("densidadeAparente_kg_L", densidadeAparente_kg_L);
  if (!(Kd_L_kg >= 0)) throw new RangeError("Kd_L_kg não pode ser negativo.");
  if (!(porosidadeTotal > 0 && porosidadeTotal <= 1)) throw new RangeError("porosidadeTotal deve estar entre 0 e 1.");
  return 1 + (densidadeAparente_kg_L * Kd_L_kg) / porosidadeTotal;
}

/** Velocidade do soluto retardado v_c = v / R. */
export function velocidadeSoluto(velocidade_m_s: number, R: number): number {
  if (!(R >= 1)) throw new RangeError("R deve ser maior ou igual a 1.");
  return velocidade_m_s / R;
}

/* ------------------------------------------------------------------ unidades */
const SEG_DIA = 86400;
export const SEGUNDOS_POR_ANO = 365.25 * SEG_DIA;

export const converterK = {
  cmPorSParaMPorS: (v: number) => v / 100,
  mPorSParaCmPorS: (v: number) => v * 100,
  mPorDiaParaMPorS: (v: number) => v / SEG_DIA,
  mPorSParaMPorDia: (v: number) => v * SEG_DIA,
  mPorSParaMPorAno: (v: number) => v * SEGUNDOS_POR_ANO,
};
