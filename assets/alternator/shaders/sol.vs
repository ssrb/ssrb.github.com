attribute vec2 aPos;
attribute float aSol;
varying float sol;

uniform mat4 mvMatrix;
uniform mat4 prMatrix;

uniform float magnitude, sign;

void main(void) {
 gl_Position = prMatrix * mvMatrix * vec4(aPos, 0.0, 1.0);
 sol = (1. + sign * aSol / magnitude) / 2.;
}