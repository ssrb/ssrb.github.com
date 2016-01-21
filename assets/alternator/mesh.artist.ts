// Copyright (c) 2016, Sebastien Sydney Robert Bigot
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies,
// either expressed or implied, of the FreeBSD Project.

///<reference path="mesh.ts" />
///<reference path="typings/gl-matrix/gl-matrix.d.ts"/>
///<reference path="typings/browserify/browserify.d.ts"/>

import msh = require('./mesh');
import Mesh = msh.Mesh;

var glmat = require('./bower_components/gl-matrix/dist/gl-matrix-min.js');

export class MeshArtist {

    public constructor(gl: WebGLRenderingContext, mesh: Mesh) {
        this.gl = gl;
        this.mesh = mesh;
        this.init();
    }

    private init() {
        var gl = this.gl;

        gl.getExtension('OES_standard_derivatives');

        // Browserify will bundle shaders and js all together for us.
        // In order to do so, the tool must find a 'require' with a string literal argument
        // to figure out what must be bundled together
        require('./shaders/mesh.vs');
        require('./shaders/mesh.fs');
        require('./shaders/sol.vs');
        require('./shaders/sol.fs');

        this.meshProgram = gl.createProgram();
        gl.attachShader(this.meshProgram, MeshArtist.getShader( gl, './shaders/mesh.vs' ));
        gl.attachShader(this.meshProgram, MeshArtist.getShader( gl, './shaders/mesh.fs' ));
        gl.linkProgram(this.meshProgram);

        this.solProgram = gl.createProgram();
        gl.attachShader(this.solProgram, MeshArtist.getShader( gl, './shaders/sol.vs' ));
        gl.attachShader(this.solProgram, MeshArtist.getShader( gl, './shaders/sol.fs' ));
        gl.linkProgram(this.solProgram);

        var vertices = this.mesh.vertices;
        var triangles = this.mesh.triangles;
        var domainIndex = this.mesh.domainIndex;
        var ntriangles = triangles.length / 3;

        this.triangles = new Array<number>(6 * ntriangles);
        this.domainIndex = new Array<number>(3 * ntriangles);
        this.baryCoordinates = new Array<number>(9 * ntriangles);
        for (var i = 0; i < ntriangles; ++i) {

            this.triangles[6 * i]       = vertices[2 * triangles[3 * i]];
            this.triangles[6 * i + 1]   = vertices[2 * triangles[3 * i] + 1];

            this.triangles[6 * i + 2]   = vertices[2 * triangles[3 * i + 1]];
            this.triangles[6 * i + 3]   = vertices[2 * triangles[3 * i + 1] + 1];

            this.triangles[6 * i + 4]   = vertices[2 * triangles[3 * i + 2]];
            this.triangles[6 * i + 5]   = vertices[2 * triangles[3 * i + 2] + 1];

            this.domainIndex[3 * i]  = this.domainIndex[3 * i + 1] = this.domainIndex[3 * i + 2] = domainIndex[i];

            this.baryCoordinates[9 * i]     = 1;
            this.baryCoordinates[9 * i + 1] = 0;
            this.baryCoordinates[9 * i + 2] = 0;

            this.baryCoordinates[9 * i + 3] = 0;
            this.baryCoordinates[9 * i + 4] = 1;
            this.baryCoordinates[9 * i + 5] = 0;

            this.baryCoordinates[9 * i + 6] = 0;
            this.baryCoordinates[9 * i + 7] = 0;
            this.baryCoordinates[9 * i + 8] = 1;
        }
    }

    public draw(prMatrix : Float32Array, mvMatrix : Float32Array) {
        var gl = this.gl;

        gl.useProgram(this.meshProgram);

        var location = gl.getAttribLocation(this.meshProgram, 'aPos')
        gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(this.triangles), gl.STATIC_DRAW);
        gl.vertexAttribPointer(location, 2, gl.FLOAT, false, 0, 0);
        gl.enableVertexAttribArray(location);

        location = gl.getAttribLocation(this.meshProgram, 'aBaryCoord')
        gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(this.baryCoordinates), gl.STATIC_DRAW);
        gl.vertexAttribPointer(location, 3, gl.FLOAT, false, 0, 0);
        gl.enableVertexAttribArray(location);

        gl.uniformMatrix4fv(gl.getUniformLocation(this.meshProgram, 'prMatrix'), false, prMatrix);
        mvMatrix = new Float32Array(mvMatrix);
        for (var i = 0; i < 12; ++i) {
            gl.uniformMatrix4fv(gl.getUniformLocation(this.meshProgram, 'mvMatrix'), false, mvMatrix);
            gl.drawArrays(gl.TRIANGLES, 0, this.triangles.length / 2);
            glmat.mat4.rotateZ(mvMatrix, mvMatrix, 2 * Math.PI / 12);
        }
    }

    public drawSol(sol: number[], magnitude: number, prMatrix : Float32Array, mvMatrix : Float32Array) {
        var gl = this.gl;
        var mesh = this.mesh;

        gl.useProgram(this.solProgram);

        var location = gl.getAttribLocation(this.solProgram, 'aPos')
        gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh.vertices), gl.STATIC_DRAW);
        gl.vertexAttribPointer(location, 2, gl.FLOAT, false, 0, 0);
        gl.enableVertexAttribArray(location);

        var solLocation = gl.getAttribLocation(this.solProgram, "aSol");
        gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(sol), gl.STATIC_DRAW);
        gl.vertexAttribPointer(solLocation, 1, gl.FLOAT, false, 0, 0);
        gl.enableVertexAttribArray( solLocation );

        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, gl.createBuffer());
        gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(mesh.triangles), gl.STATIC_DRAW);

        gl.uniform1f(gl.getUniformLocation(this.solProgram, 'magnitude'), magnitude);

        gl.uniformMatrix4fv(gl.getUniformLocation(this.solProgram, 'prMatrix'), false, prMatrix);
        mvMatrix = new Float32Array(mvMatrix);
        for (var i = 0; i < 12; ++i) {
            gl.uniformMatrix4fv(gl.getUniformLocation(this.solProgram, 'mvMatrix'), false, mvMatrix);
            gl.uniform1f(gl.getUniformLocation(this.solProgram, 'sign'), i & 1 ? -1 : 1);
            gl.drawElements(gl.TRIANGLES, mesh.triangles.length, gl.UNSIGNED_SHORT ,0);
            glmat.mat4.rotateZ(mvMatrix, mvMatrix, 2 * Math.PI / 12);
        }
    }

    private static getShader(gl : WebGLRenderingContext, path : string) : WebGLShader {

        var shader : WebGLShader;

        var ext = path.substring(path.lastIndexOf(".") + 1);

        if ( ext == 'fs' )
            shader = gl.createShader ( gl.FRAGMENT_SHADER );
        else if ( ext == 'vs' )
            shader = gl.createShader(gl.VERTEX_SHADER);
        else return null;

        var glsl = require(path);

        gl.shaderSource(shader, glsl());
        gl.compileShader(shader);
        if (gl.getShaderParameter(shader, gl.COMPILE_STATUS) == 0)
            alert(path + "\n" + gl.getShaderInfoLog(shader));

        return shader;
    }

    private gl: WebGLRenderingContext;
    private mesh: Mesh;
    private meshProgram : WebGLProgram;
    private solProgram : WebGLProgram;
    private triangles : number[];
    private domainIndex : number[];
    private baryCoordinates : number[];
}