import { describe, expect, it } from "vitest";
import {
  E1,
  ajustarRecuperacao,
  ajustarTheis,
  bouwerRice,
  coeficienteParticao,
  coeficientesBouwerRice,
  converterK,
  cooperJacob,
  fatorRetardamento,
  funcaoPocoHantush,
  gradienteHidraulico,
  hvorslev,
  lnReRw,
  rebaixamentoHantush,
  rebaixamentoTheis,
  recuperacaoTheis,
  regressaoLinear,
  tempoBasicoHvorslev,
  tempoTransito,
  velocidadeLinear,
  velocidadeSoluto,
} from "./index";

const rel = (a: number, b: number) => Math.abs(a - b) / Math.abs(b);

describe("integral exponencial E1 = W(u) de Theis", () => {
  // Valores de referência: Abramowitz & Stegun, Tabela 5.1
  it.each([
    [0.01, 4.0379295765],
    [0.1, 1.8229239584],
    [0.5, 0.5597735948],
    [1, 0.2193839344],
    [2, 0.04890051071],
    [5, 0.001148295591],
  ])("E1(%f) = %f", (x, ref) => {
    expect(rel(E1(x), ref)).toBeLessThan(1e-8);
  });
  it("é contínua na troca de método (x = 1)", () => {
    expect(rel(E1(1 - 1e-9), E1(1 + 1e-9))).toBeLessThan(1e-6);
  });
});

describe("Hantush-Jacob", () => {
  it("β = 0 coincide com Theis", () => {
    expect(funcaoPocoHantush(0.01, 0)).toBe(E1(0.01));
  });
  it("para u muito pequeno tende a 2·K0(β) (regime permanente de De Glee)", () => {
    // K0(0,1) = 2,4270690 → 2·K0 = 4,8541380 (Abramowitz & Stegun, Tabela 9.8)
    expect(rel(funcaoPocoHantush(1e-7, 0.1), 4.854138)).toBeLessThan(2e-3);
  });
  it("a drenança reduz o rebaixamento em relação a Theis", () => {
    const base = { Q_m3_s: 0.01, T_m2_s: 1e-3, S: 1e-4, r_m: 20, t_s: 3600 };
    expect(rebaixamentoHantush({ ...base, B_m: 200 })).toBeLessThan(rebaixamentoTheis(base));
  });
});

describe("Theis, Cooper-Jacob e recuperação", () => {
  const Q = 0.01; // m³/s (36 m³/h)
  const T = 5e-3; // m²/s
  const S = 2e-4;
  const r = 30; // m
  const tempos = [60, 120, 300, 600, 1200, 1800, 3600, 7200, 14400, 28800, 43200, 86400];
  const dados = tempos.map((t) => ({ t_s: t, s_m: rebaixamentoTheis({ Q_m3_s: Q, T_m2_s: T, S, r_m: r, t_s: t }) }));

  it("o ajuste de Theis recupera T e S de dados sintéticos", () => {
    const aj = ajustarTheis(dados, { Q_m3_s: Q, r_m: r });
    expect(rel(aj.T_m2_s, T)).toBeLessThan(1e-4);
    expect(rel(aj.S, S)).toBeLessThan(1e-3);
    expect(aj.rmse_m).toBeLessThan(1e-6);
  });

  it("Cooper-Jacob nos tempos longos (u < 0,01) fica perto de Theis", () => {
    const cj = cooperJacob(dados, { Q_m3_s: Q, r_m: r, tMinimo_s: 7200 });
    expect(rel(cj.T_m2_s, T)).toBeLessThan(0.02);
    expect(rel(cj.S, S)).toBeLessThan(0.05);
    expect(cj.avisos).toHaveLength(0);
  });

  it("Cooper-Jacob avisa quando usa tempos curtos (u > 0,01)", () => {
    const cj = cooperJacob(dados, { Q_m3_s: Q, r_m: r });
    expect(cj.avisos.length).toBeGreaterThan(0);
  });

  it("recuperação de Theis: reta s′ × log(t/t′) devolve T", () => {
    const tp = 86400; // bombeou 1 dia
    const rec = [60, 300, 900, 3600, 7200, 14400].map((tl) => {
      const t = tp + tl;
      const sres =
        rebaixamentoTheis({ Q_m3_s: Q, T_m2_s: T, S, r_m: r, t_s: t }) - rebaixamentoTheis({ Q_m3_s: Q, T_m2_s: T, S, r_m: r, t_s: tl });
      return { t_s: t, tLinha_s: tl, sResidual_m: sres };
    });
    const res = recuperacaoTheis(rec.slice(3), Q); // t′ grande: aproximação válida
    expect(rel(res.T_m2_s, T)).toBeLessThan(0.03);
  });
});

describe("slug test", () => {
  const T0 = 120; // s
  const leituras = [0, 10, 20, 40, 60, 90, 120, 180, 240, 300].map((t) => ({ t_s: t, h_m: 0.8 * Math.exp(-t / T0) }));

  it("encontra o tempo básico T0 de Hvorslev", () => {
    const aj = ajustarRecuperacao(leituras, [0.05, 1]);
    expect(rel(tempoBasicoHvorslev(aj), T0)).toBeLessThan(1e-9);
    expect(aj.reta.r2).toBeGreaterThan(0.999999);
  });

  it("calcula K de Hvorslev pela fórmula de Le/R > 8", () => {
    const rc = 0.025, R = 0.05, Le = 1.5;
    const { K_m_s, avisos } = hvorslev({ rc_m: rc, R_m: R, Le_m: Le, T0_s: T0 });
    expect(rel(K_m_s, (rc * rc * Math.log(Le / R)) / (2 * Le * T0))).toBeLessThan(1e-12);
    expect(K_m_s).toBeGreaterThan(1e-6);
    expect(K_m_s).toBeLessThan(1e-5);
    expect(avisos).toHaveLength(0);
  });

  it("avisa quando Le/R ≤ 8", () => {
    expect(hvorslev({ rc_m: 0.025, R_m: 0.1, Le_m: 0.5, T0_s: 60 }).avisos.length).toBe(1);
  });

  it("coeficientes A, B, C de Bouwer & Rice ficam nas curvas originais", () => {
    // Leitura aproximada das curvas de Bouwer & Rice (1976): Le/rw = 10 → A ≈ 1,8, B ≈ 0,3, C ≈ 1,2;
    // Le/rw = 100 → A ≈ 4,3, B ≈ 0,7, C ≈ 4,2
    const c10 = coeficientesBouwerRice(10);
    const c100 = coeficientesBouwerRice(100);
    expect(c10.A).toBeCloseTo(1.82, 1);
    expect(c10.B).toBeCloseTo(0.29, 1);
    expect(c10.C).toBeCloseTo(1.19, 1);
    expect(c100.A).toBeCloseTo(4.29, 1);
    expect(c100.B).toBeCloseTo(0.72, 1);
    expect(c100.C).toBeCloseTo(4.25, 1);
  });

  it("ln(Re/rw) cai na faixa usual (≈ 1,5 a 4) e K sai positivo", () => {
    const geo = { rw_m: 0.075, Le_m: 3, Lw_m: 4, H_m: 15 };
    const ln = lnReRw(geo);
    expect(ln.valor).toBeGreaterThan(1.5);
    expect(ln.valor).toBeLessThan(4);
    const aj = ajustarRecuperacao(leituras, [0.05, 1]);
    const br = bouwerRice({ rc_m: 0.025, ...geo, inclinacao_por_s: aj.reta.b });
    expect(br.K_m_s).toBeGreaterThan(0);
    // poço completo (Lw = H) usa o coeficiente C
    expect(lnReRw({ rw_m: 0.075, Le_m: 3, Lw_m: 15, H_m: 15 }).valor).toBeGreaterThan(0);
  });
});

describe("Darcy, retardamento e tempo de trânsito", () => {
  it("velocidade linear média e tempo até o receptor", () => {
    const i = gradienteHidraulico(101.5, 100, 300); // 0,005
    expect(i).toBeCloseTo(0.005, 12);
    const v = velocidadeLinear(1e-4, i, 0.25); // 2e-6 m/s
    expect(v).toBeCloseTo(2e-6, 15);
    expect(converterK.mPorSParaMPorAno(v)).toBeCloseTo(63.1152, 3);
    expect(tempoTransito(150, v) / 86400).toBeCloseTo(868.0556, 3); // dias
  });

  it("fator de retardamento com Koc·foc (benzeno, exemplo)", () => {
    const Kd = coeficienteParticao(59, 0.002); // 0,118 L/kg
    const R = fatorRetardamento(1.6, Kd, 0.3);
    expect(R).toBeCloseTo(1.629333, 5);
    expect(velocidadeSoluto(2e-6, R)).toBeCloseTo(2e-6 / R, 15);
  });

  it("recusa entradas fisicamente impossíveis", () => {
    expect(() => velocidadeLinear(1e-4, 0.01, 0)).toThrow();
    expect(() => coeficienteParticao(59, 2)).toThrow();
    expect(() => fatorRetardamento(1.6, -1, 0.3)).toThrow();
  });

  it("regressão linear básica", () => {
    const r = regressaoLinear([0, 1, 2, 3], [1, 3, 5, 7]);
    expect(r.a).toBeCloseTo(1, 12);
    expect(r.b).toBeCloseTo(2, 12);
    expect(r.r2).toBeCloseTo(1, 12);
  });
});
