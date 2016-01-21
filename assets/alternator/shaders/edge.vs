attribute vec2 aPos;
attribute float aBoundaryIndex;
uniform mat4 mvMatrix, prMatrix;
varying float vBoundaryIndex;
void main(void) {
    gl_Position = prMatrix * mvMatrix * vec4(aPos, 0.0, 1.0);
    vBoundaryIndex = aBoundaryIndex;
}