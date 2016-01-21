precision highp float;
varying float sol;

vec3 heatcolor() {

  float h = 4. * sol,
  s = 1.,
  v = 0.8,
  f = h - floor(h),
  p = v * (1. - s),
  q = v * (1. - s * f),
  t = v * (1. - s * (1. - f));

  if (h <= 1.) {
    return vec3(v, t, p);
  }

  if (h <= 2.) {
    return vec3(q, v, p);
  }

  if (h <= 3.) {
    return vec3(p, v, t);
  }

  return vec3(p, q, v);
}

void main(void) {
   gl_FragColor = vec4(heatcolor(), 1.);
}