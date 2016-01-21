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

///<reference path="mesh.ts"/>
///<reference path="my_typings/numericjs/numericjs.d.ts"/>
///<reference path="node_modules/typescript-collections/collections.d.ts"/>

require('./bower_components/numericjs/lib/numeric-1.2.6.min.js');

import msh = require('./mesh');
import Mesh = msh.Mesh;

import rct = require('./rct');
import ReverseConnectivityTable = rct.ReverseConnectivityTable;

export enum BoundaryType {
    Interior = 0,
    Sliding = 1,
    Up = 2,
    Down = 3,
    Centre = 4,
    Outside = 5
}

export enum DomainType {
    Air = 0,
    RotorIron = 1,
    StatorIron = 2,
    RotorCopper = 3,
    SupplyCoilPositive = 19,
    InductorCoilPositive = 20,
    SupplyCoilNegative = 21,
    InductorCoilNegative = 22
}

export class Domain {

    public constructor(mesh : Mesh) {

        this.mesh = mesh;

        this.rct = new ReverseConnectivityTable(mesh);
        this.indexPhases();
        this.indexBoundaries();
        this.precomputeTrisArea();

    };

    private indexBoundaries() {
        var mesh = this.mesh;
        this.up = new collections.Set<number>();
        this.down = new collections.Set<number>();
        this.sliding = new collections.Set<number>();
        for (var ei = 0; ei < mesh.edges.length / 2; ++ei) {
            var vi = mesh.edges[2 * ei], vj = mesh.edges[2 * ei + 1];
            switch (mesh.boundaryIndex[ei]) {
                case BoundaryType.Up:
                    this.up.add(vi);
                    this.up.add(vj);
                    break;
                case BoundaryType.Down:
                    this.down.add(vi);
                    this.down.add(vj);
                    break;
                case BoundaryType.Sliding:
                    this.sliding.add(vi);
                    this.sliding.add(vj);
                    break;
                default:
                    break;
            }
        }
    };

    private indexPhases() {
        var tris = this.mesh.triangles;
        var ntris = tris.length / 3;
        this.phases = new collections.Set<number>();
        for (var ti = 0; ti < ntris; ++ti) {
            switch (this.mesh.domainIndex[ti]) {
                case DomainType.SupplyCoilPositive:
                case DomainType.InductorCoilPositive:
                case DomainType.SupplyCoilNegative:
                case DomainType.InductorCoilNegative:
                    this.phases.add(tris[3 * ti]);
                    this.phases.add(tris[3 * ti + 1]);
                    this.phases.add(tris[3 * ti + 2]);
                    break;
                default:
                    break;
            }
        }
    };

    private precomputeTrisArea() {
        var mesh = this.mesh;
        var verts = mesh.vertices;
        var tris = mesh.triangles;
        var ntris = tris.length / 3;
        this.q = numeric.rep([ntris, 3], []);
        this.area = numeric.rep([ntris], 0);
        for (var ti = 0; ti < ntris; ++ti) {
            for (var si = 0; si < 3; ++si) {
                var xi = 3 * ti + si, xj = 3 * ti + ((si + 1) % 3);
                // Mesh is in mm but we use MKSA
                this.q[ti][si] = [
                    (verts[2 * tris[xi]] - verts[2 * tris[xj]]) / 1000,
                    (verts[2 * tris[xi] + 1] - verts[2 * tris[xj] + 1]) / 1000
                ];
            }
            this.area[ti] = 0.5 * numeric.det([this.q[ti][0], this.q[ti][1]]);
        }
    };

    public applyAntiPeriodicBoundaryConditions(rotation: number) : void {

        var mesh = this.mesh;
        var up = this.up;
        var down = this.down;

        var verts = mesh.vertices;
        var nverts = verts.length / 2;

        var v2dof = new Array<number>(nverts);
        var coeff = numeric.rep([nverts], 1);

        var ndof = 0;
        for (var vi = 0; vi < nverts; ++vi) {
            if (!this.up.contains(vi) && !this.down.contains(vi))
            {
                v2dof[vi] = ndof++;
            }
        }

        function hypot (x: number, y: number) : number {
            return Math.sqrt(x * x + y * y);
        }

        var rotorUpTouchStator = rotation <= 0;

        var centre = -1;
        down.forEach( function(vd) {
            v2dof[vd] = ndof;
            var downLen = hypot(verts[2 * vd], verts[2 * vd + 1]);
            up.forEach(function (vu : number) {
                if (vd != vu) {
                    var upLen = hypot(verts[2 * vu], verts[2 * vu + 1]);
                    if (Math.abs(downLen - upLen) < Domain.kEpsilon) {
                        v2dof[vu] = ndof;
                        if (rotorUpTouchStator) {
                            coeff[vd] = -coeff[vu];
                        } else {
                            coeff[vu] = -coeff[vd];
                        }
                        return false;
                    }
                } else {
                    centre = vu;
                    return false;
                }
            });
            ++ndof;
        });

        this.ndof = ndof;
        this.v2dof = v2dof;
        this.coeff = coeff;
        this.centre = centre;
    };

    public static joinSlidingDomains(rotor: Domain, stator : Domain, rotation: number) : number {

        var ndof = 0;

        function mapLocalToGlobal(domain : Domain) : number[] {
            var nverts = domain.mesh.vertices.length / 2;
            var l2g = numeric.rep([domain.ndof], -1);
            for (var vi = 0; vi < nverts; ++vi) {
                var di = domain.v2dof[vi];
                if (!domain.sliding.contains(vi) && l2g[di] == -1) {
                    l2g[di] = ndof++;
                }
            }
            return l2g;
        }

        function applyMapping(data : number[], mapping : number[]) : void {
            for (var vi = 0; vi < data.length; ++vi) {
                data[vi] = mapping[data[vi]];
            }
        }

        var sl2g = mapLocalToGlobal(stator);
        var rl2g = mapLocalToGlobal(rotor);

        var rverts = rotor.mesh.vertices, sverts = stator.mesh.vertices;
        rotor.sliding.forEach(function(vr) {
            var dr = rotor.v2dof[vr];
            if (rl2g[dr] == -1) {
                rl2g[dr] = ndof;
                var rtheta = Math.atan2(rverts[2 * vr + 1], rverts[2 * vr]);
                rtheta += rotation;
                var sign = 1;
                if (rtheta > Math.PI / 12) {
                    rtheta -= Math.PI / 6;
                    sign *= -1;
                }
                if (rtheta < -Math.PI / 12) {
                    rtheta += Math.PI / 6;
                    sign *= -1;
                }
                stator.sliding.forEach(function (vs) {
                    var ds = stator.v2dof[vs];
                    var stheta = Math.atan2(sverts[2 * vs + 1], sverts[2 * vs]);
                    if (Math.abs(rtheta - stheta) < Domain.kEpsilon && sl2g[ds] == -1) {
                        sl2g[ds] = ndof;
                        rotor.coeff[vr] = sign * stator.coeff[vs];
                        return false;
                    }
                });
                ++ndof;
            }
        });

        applyMapping(rotor.v2dof, rl2g);
        applyMapping(stator.v2dof, sl2g);

        return ndof;
    };

    static kEpsilon = 10e-6;

    mesh : Mesh;
    ndof : number;
    v2dof : number[];
    coeff : number[];
    up : collections.Set<number>;
    down : collections.Set<number>;
    sliding : collections.Set<number>;
    phases : collections.Set<number>;
    centre : number;
    rct : ReverseConnectivityTable;
    area: number[];
    q: number[][][];
}

