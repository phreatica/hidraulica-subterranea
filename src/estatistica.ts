/** Regressão linear por mínimos quadrados y = a + b·x, com R². */
export interface Reta {
  a: number; // intercepto
  b: number; // inclinação
  r2: number;
  n: number;
}

export function regressaoLinear(x: number[], y: number[]): Reta {
  if (x.length !== y.length) throw new RangeError("x e y precisam ter o mesmo tamanho.");
  const n = x.length;
  if (n < 2) throw new RangeError("São necessários pelo menos 2 pontos.");
  let sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sx += x[i]; sy += y[i]; sxx += x[i] * x[i]; sxy += x[i] * y[i]; syy += y[i] * y[i];
  }
  const dx = n * sxx - sx * sx;
  if (dx === 0) throw new RangeError("Todos os x são iguais: não há reta.");
  const b = (n * sxy - sx * sy) / dx;
  const a = (sy - b * sx) / n;
  const dy = n * syy - sy * sy;
  const r = dy === 0 ? 1 : (n * sxy - sx * sy) / Math.sqrt(dx * dy);
  return { a, b, r2: r * r, n };
}
