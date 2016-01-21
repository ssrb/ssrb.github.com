declare class LU
{
    L : Array<number>;
    U : Array<number>;
    P : Array<number>;
    Pinv : Array<number>;
}

declare module numeric {
    export function ccsSparse(A : number[][]): [number[], number[], number[]];
    export function ccsLUP(A : [number[], number[], number[]], threshold : number) : LU;
    // return {L:L, U:U, P:P, Pinv:Pinv};

    export function ccsLUPSolve(LUP : LU, B : number[]) : number[];
    // return b;
    // return [Xi,Xj,Xv];

    export function det(x : [number[], number[]]) : number;

    export function dot(x : any[], y : any[]) : any;

    export function rep<T>(dimensions : number[], init : T) : any;
}
