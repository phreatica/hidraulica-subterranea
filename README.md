# hidraulica-subterranea

[![licença MIT](https://img.shields.io/badge/licen%C3%A7a-MIT-blue)](LICENSE)

Hidráulica de águas subterrâneas em TypeScript, sem dependências: **slug test**, **ensaio de bombeamento**, **Darcy**, **retardamento** e **tempo de trânsito**. Unidades SI, fórmulas citadas no código e testes com valores de referência.

É o motor comum dos produtos do ecossistema [Phreatica](https://phreatica.com) e [Estrato Perfis](https://estratoperfis.com.br): o mesmo cálculo, com o mesmo resultado, nos dois.

## O que tem

| Grupo | Funções |
|---|---|
| **Slug test** | `ajustarRecuperacao` (reta ln(h/h₀) × t), `tempoBasicoHvorslev`, `hvorslev`, `bouwerRice`, `lnReRw`, `coeficientesBouwerRice` (A, B, C pela aproximação de Van Rooy, 1988) |
| **Bombeamento** | `E1` / `funcaoPocoTheis` W(u), `rebaixamentoTheis`, `ajustarTheis` (mínimos quadrados), `cooperJacob` (com aviso quando u > 0,01), `recuperacaoTheis`, `funcaoPocoHantush` W(u, r/B), `rebaixamentoHantush` |
| **Transporte** | `gradienteHidraulico`, `fluxoDarcy`, `velocidadeLinear`, `tempoTransito`, `coeficienteParticao` (K_oc·f_oc), `fatorRetardamento`, `velocidadeSoluto`, `converterK` (cm/s, m/dia, m/ano) |

## Exemplo

```ts
import { cooperJacob, velocidadeLinear, fatorRetardamento, coeficienteParticao, converterK } from "hidraulica-subterranea";

// Cooper-Jacob: só os tempos em que u < 0,01
const ensaio = cooperJacob(leituras, { Q_m3_s: 0.01, r_m: 30, tMinimo_s: 7200 });
console.log(ensaio.T_m2_s, ensaio.S, ensaio.avisos);

// Velocidade da água e do benzeno até um receptor
const v = velocidadeLinear(1e-4, 0.005, 0.25);                 // m/s
const R = fatorRetardamento(1.6, coeficienteParticao(59, 0.002), 0.3);
console.log(converterK.mPorSParaMPorAno(v), converterK.mPorSParaMPorAno(v / R)); // m/ano
```

## Princípio

A biblioteca calcula; **a interpretação é do profissional responsável**. As funções devolvem avisos objetivos (por exemplo, "u = 0,03 no primeiro ponto: acima de 0,01") em vez de decidir pelo usuário.

## Referências
- Hvorslev, M. J. (1951). *Time lag and soil permeability in ground-water observations*. USACE Bulletin 36.
- Bouwer, H.; Rice, R. C. (1976). *Water Resources Research*, 12(3), 423–428. Bouwer, H. (1989). *Ground Water*, 27(3), 304–309.
- Van Rooy, D. (1988). Aproximação polinomial das curvas A, B e C de Bouwer & Rice.
- Theis, C. V. (1935). *Transactions AGU*, 16, 519–524. Cooper, H. H.; Jacob, C. E. (1946). *Transactions AGU*, 27, 526–534.
- Hantush, M. S.; Jacob, C. E. (1955). *Transactions AGU*, 36, 95–100.
- Kruseman, G. P.; De Ridder, N. A. (1990). *Analysis and Evaluation of Pumping Test Data*. ILRI, Publ. 47.
- Abramowitz, M.; Stegun, I. A. (1964). *Handbook of Mathematical Functions* (valores de E₁ e K₀ usados nos testes).

## English
Groundwater hydraulics toolkit in TypeScript: slug tests (Hvorslev, Bouwer & Rice), pumping tests (Theis, Cooper-Jacob, Theis recovery, Hantush-Jacob), Darcy velocity, retardation factor and travel time. SI units, no dependencies, tested against published reference values. API names are in Portuguese.

## Licença
MIT © 2026 [Phreatica](https://phreatica.com).
