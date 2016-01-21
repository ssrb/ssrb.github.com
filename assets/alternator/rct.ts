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

///<reference path="my_typings/numericjs/numericjs.d.ts"/>
///<reference path="mesh.ts"/>

require('./bower_components/numericjs/lib/numeric-1.2.6.min.js');

import msh = require('./mesh');
import Mesh = msh.Mesh;

interface ILoopFunction {
    (ti : number): void;
}

export class ReverseConnectivityTable {

    public constructor(mesh : Mesh) {
        var nverts = mesh.vertices.length / 2;
        this.head = numeric.rep([nverts], -1);
        this.next = numeric.rep([mesh.triangles.length], -1);
        for (var p = 0; p < mesh.triangles.length; ++p) {
            var vi = mesh.triangles[p];
            this.next[p] = this.head[vi];
            this.head[vi] = p; // (p / 3) = ti, is the triangle number, the new head of the list of triangles for vertex vi;
        }
    }

    public forEach(vi : number, callback: ILoopFunction) {
        var p = this.head[vi];
        while (p != -1) {
            callback(Math.floor(p / 3));
            p = this.next[p];
        }
    }

    private head : number[];
    private next : number[];
}