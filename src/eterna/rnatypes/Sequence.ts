import EPars, {RNABase} from 'eterna/EPars';
import SecStruct from './SecStruct';

export default class Sequence {
    constructor(seq: string) {
        this._sequence = seq;
    }

    public hasCut(from: number, to: number): boolean {
        return this.sequence.slice(from, to + 1).some(
            (c) => c === RNABase.CUT
        );
    }

    public static fromSequence(seq: number[]) {
        return new Sequence(EPars.sequenceToString(seq));
    }

    public getColoredSequence(): string {
        let res = '';
        for (const c of this._sequence) {
            res += EPars.getColoredLetter(c);
        }
        return res;
    }

    public getExpColoredSequence(expData: number[]): string {
        // AMW TODO: how could this be?
        if (expData == null) {
            return this._sequence;
        }

        const offset: number = expData[0];
        const maxmax: number = Math.max(...expData.slice(1));
        const minmin: number = Math.min(...expData.slice(1));

        const avg: number = (maxmax + minmin) / 2.0;

        let res = '';
        for (let ii = 0; ii < this._sequence.length; ii++) {
            if (ii < offset - 1 || ii >= expData.length) {
                res += this._sequence[ii];
            } else if (expData[ii] < avg) {
                res += `<FONT COLOR='#7777FF'>${this._sequence[ii]}</FONT>`;
            } else {
                res += `<FONT COLOR='#FF7777'>${this._sequence[ii]}</FONT>`;
            }
        }

        return res;
    }

    public countConsecutive(letter: number, locks: boolean[] | null = null): number {
        const sequence = EPars.stringToSequence(this._sequence);

        let maxConsecutive = 0;

        let ii = 0;
        let startIndex = -1;
        for (ii = 0; ii < sequence.length; ii++) {
            if (sequence[ii] === letter) {
                if (startIndex < 0) {
                    startIndex = ii;
                }
            } else if (startIndex >= 0) {
                if (maxConsecutive < ii - startIndex) {
                    if (locks == null) {
                        maxConsecutive = ii - startIndex;
                    } else {
                        let allLocked = true;
                        let jj: number;
                        for (jj = startIndex; jj < ii; jj++) {
                            allLocked = allLocked && locks[jj];
                        }
                        if (allLocked === false) {
                            maxConsecutive = ii - startIndex;
                        }
                    }
                }
                startIndex = -1;
            }
        }

        if (startIndex >= 0) {
            if (maxConsecutive < ii - startIndex) {
                maxConsecutive = ii - startIndex;
            }
        }

        return maxConsecutive;
    }

    public getRestrictedConsecutive(
        letter: number, maxAllowed: number, locks: boolean[] | null = null
    ): number[] {
        const sequence = EPars.stringToSequence(this._sequence);

        const restricted: number[] = [];

        let ii = 0;
        let startIndex = -1;

        if (maxAllowed <= 0) {
            return restricted;
        }

        for (ii = 0; ii < sequence.length; ii++) {
            if (sequence[ii] === letter) {
                if (startIndex < 0) {
                    startIndex = ii;
                }
            } else if (startIndex >= 0) {
                if (maxAllowed < ii - startIndex) {
                    if (locks == null) {
                        restricted.push(startIndex);
                        restricted.push(ii - 1);
                    } else {
                        let allLocked = true;
                        let jj: number;
                        for (jj = startIndex; jj < ii; jj++) {
                            allLocked = allLocked && locks[jj];
                        }
                        if (allLocked === false) {
                            restricted.push(startIndex);
                            restricted.push(ii - 1);
                        }
                    }
                }
                startIndex = -1;
            }
        }

        // gotta check if we found a startIndex without an end...
        if (startIndex >= 0) {
            if (maxAllowed < ii - startIndex) {
                restricted.push(startIndex);
                restricted.push(ii - 1);
            }
        }

        return restricted;
    }

    public getSequenceRepetition(n: number): number {
        const dict: Set<string> = new Set<string>();
        let numRepeats = 0;

        for (let ii = 0; ii < this._sequence.length - n; ii++) {
            const substr: string = this._sequence.substr(ii, n);
            if (dict.has(substr)) {
                numRepeats++;
            } else {
                dict.add(substr);
            }
        }

        return numRepeats++;
    }

    public numGUPairs(pairs: SecStruct): number {
        const sequence = EPars.stringToSequence(this._sequence);
        let ret = 0;

        for (let ii = 0; ii < pairs.length; ii++) {
            if (pairs.pairs[ii] > ii) {
                if (sequence[ii] === RNABase.GUANINE && sequence[pairs.pairs[ii]] === RNABase.URACIL) {
                    ret++;
                }
                if (sequence[ii] === RNABase.URACIL && sequence[pairs.pairs[ii]] === RNABase.GUANINE) {
                    ret++;
                }
            }
        }

        return ret;
    }

    public numGCPairs(pairs: SecStruct): number {
        const sequence = EPars.stringToSequence(this._sequence);
        let ret = 0;

        for (let ii = 0; ii < pairs.length; ii++) {
            if (pairs.pairs[ii] > ii) {
                if (sequence[ii] === RNABase.GUANINE && sequence[pairs.pairs[ii]] === RNABase.CYTOSINE) {
                    ret++;
                }
                if (sequence[ii] === RNABase.CYTOSINE && sequence[pairs.pairs[ii]] === RNABase.GUANINE) {
                    ret++;
                }
            }
        }

        return ret;
    }

    public numUAPairs(pairs: SecStruct): number {
        const sequence = EPars.stringToSequence(this._sequence);
        let ret = 0;

        for (let ii = 0; ii < pairs.length; ii++) {
            if (pairs.pairs[ii] > ii) {
                if (sequence[ii] === RNABase.ADENINE && sequence[pairs.pairs[ii]] === RNABase.URACIL) {
                    ret++;
                }
                if (sequence[ii] === RNABase.URACIL && sequence[pairs.pairs[ii]] === RNABase.ADENINE) {
                    ret++;
                }
            }
        }

        return ret;
    }

    public get sequence(): number[] {
        return EPars.stringToSequence(this._sequence);
    }

    public set sequence(sequence: number[]) {
        this._sequence = EPars.sequenceToString(sequence);
    }

    public get sequenceString(): string {
        return this._sequence;
    }

    public get length(): number {
        return this._sequence.length;
    }

    /**
     * Returns a new sequence using indices that obey slice/substr notation, i.e.,
     * start to end-1 are cut out
     * @param start The first index to begin slicing
     * @param end The one-past-the end index; with start === 0, doubles as len
     */
    public slice(start: number, end: number = -1): Sequence {
        if (end === -1) {
            return new Sequence(this._sequence.substr(start));
        } else {
            return new Sequence(this._sequence.substr(start, end));
        }
    }

    private _sequence: string;
}
