#extension GL_OES_standard_derivatives : enable
precision highp float;
varying vec3 vBaryCoord;

float edgeFactor(){
  vec3 d = fwidth(vBaryCoord);
  vec3 a3 = smoothstep(vec3(0.0), d * 0.1, vBaryCoord);
  return 1.0 - min(min(a3.x, a3.y), a3.z);
}

void main(void) {
  gl_FragColor = vec4(1.0, 1.0, 1.0, edgeFactor() * 0.5);
}