uniform mat4 mvMatrix, prMatrix;
attribute vec2 aPos;
attribute vec3 aBaryCoord;
varying vec3 vBaryCoord;
void main(void) {
   gl_Position = prMatrix * mvMatrix * vec4(aPos, 0.0, 1.0);
   vBaryCoord = aBaryCoord;
}