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
export class Mesh {

    constructor(vertices: number[],
                triangles: number[],
                domainIndex: number[],
                edges: number[],
                boundaryIndex: number[]) {
        this.vertices = vertices;
        this.triangles = triangles;
        this.domainIndex = domainIndex;
        this.edges = edges;
        this.boundaryIndex = boundaryIndex;
    };

    public static load(req:XMLHttpRequest):Mesh {
        if (req.readyState == XMLHttpRequest.DONE) {

            var lines = req.responseText.split("\n");
            var array = lines[0].match(/^\s*(\S+)\s+(\S+)\s+(\S+)$/);
            var nverts = parseInt(array[1]);
            var ntris = parseInt(array[2]);
            var nedges = parseInt(array[3]);

            var li = 1;
            var vertices = new Array<number>(2 * nverts);
            for (var vi = 0; vi < nverts; ++vi) {
                array = lines[li++].match(/^\s*(\S+)\s*(\S+)\s+(\S+)$/);
                vertices[2 * vi] = parseFloat(array[1]) / 300;
                vertices[2 * vi + 1] = parseFloat(array[2]) / 300;
            }

            var triangles = new Array<number>(3 * ntris);
            var domainIndex = new Array<number>(ntris);
            for (var ti = 0; ti < ntris; ++ti) {
                array = lines[li++].match(/^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)$/);

                var i = parseInt(array[1]) - 1, j = parseInt(array[2]) - 1, k = parseInt(array[3]) - 1, d = parseInt(array[4]);

                triangles[3 * ti] = i;
                triangles[3 * ti + 1] = j;
                triangles[3 * ti + 2] = k;

                domainIndex[ti] = d;
            }

            var edges = new Array<number>(2 * nedges);
            var boundaryIndex = new Array<number>(nedges);
            for (var ei = 0; ei < nedges; ++ei) {
                array = lines[li++].match(/^\s*(\S+)\s+(\S+)\s+(\S+)$/);

                var i = parseInt(array[1]) - 1, j = parseInt(array[2]) - 1;

                edges[2 * ei] = i;
                edges[2 * ei + 1] = j;

                boundaryIndex[ei] = parseInt(array[3]);
            }

            return new Mesh(vertices, triangles, domainIndex, edges, boundaryIndex);
        }
        return null;
    };

    vertices:number[];
    triangles:number[];
    domainIndex:number[];
    edges:number[];
    boundaryIndex:number[];
}


