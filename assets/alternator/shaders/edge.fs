precision highp float;
varying float vBoundaryIndex;
void main(void) {
    gl_FragColor = vBoundaryIndex < 1.5 ? vec4(1.0,0.0,0.0,1.0) : vec4(0.0,0.0,1.0,1.0);
}