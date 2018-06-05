import {EPars} from "./EPars";
import {Folder} from "./folding/Folder";
import {Plot, PlotType} from "./Plot";
import {Pose2D} from "./pose2D/Pose2D";

export enum UndoBlockParam {
    PARAM_GU = 0,
    PARAM_GC = 1,
    PARAM_AU = 2,
    PARAM_MFE = 3,
    PARAM_REPETITION = 4,
    PARAM_STACK = 5,
    PARAM_FE = 6,
    PARAM_DOTPLOT = 7,
    PARAM_DOTPLOT_BITMAP = 8,
    PARAM_MELTPLOT_BITMAP = 9,
    PARAM_PROB_SCORE = 10,
    PARAM_MELTING_POINT = 11,
    PARAM_PAIR_SCORE = 12,
    PARAM_NNFE_ARRAY = 13,
    PARAM_MAX = 14,
}

export class UndoBlock {
    public constructor(seq: number[]) {
        this._sequence = seq.slice();
    }

    public toJson(): any {
        return {
            _sequence: this._sequence,
            _pairs_array: this._pairs_array,
            _params_array: this._params_array,
            _stable: this._stable,
            _target_oligo: this._target_oligo,
            _target_oligos: this._target_oligos,
            _oligo_order: this._oligo_order,
            _oligos_paired: this._oligos_paired,
            _target_pairs: this._target_pairs,
            _target_oligo_order: this._target_oligo_order,
            _puzzle_locks: this._puzzle_locks,
            _forced_struct: this._forced_struct,
            _target_conditions: this._target_conditions
        };
    }

    public fromJson(json: any): void {
        try {
            this._sequence = json['_sequence'];
            this._pairs_array = json['_pairs_array'];
            this._params_array = json['_params_array'];
            this._stable = json['_stable'];
            this._target_oligo = json['_target_oligo'];
            this._target_oligos = json['_target_oligos'];
            this._oligo_order = json['_oligo_order'];
            this._oligos_paired = json['_oligos_paired'];
            this._target_pairs = json['_target_pairs'];
            this._target_oligo_order = json['_target_oligo_order'];
            this._puzzle_locks = json['_puzzle_locks'];
            this._forced_struct = json['_forced_struct'];
            this._target_conditions = json['_target_conditions'];
        } catch (e) {
            throw new Error(`Error parsing UndoBlock JSON: ${e}`);
        }
    }

    public get_target_oligos(): any[] {
        return this._target_oligos;
    }

    public get_target_oligo(): any[] {
        return this._target_oligo;
    }

    public get_oligo_mode(): number {
        let tc: any = this.get_target_conditions();
        if (tc == null) return 0;
        return tc['fold_mode'] == null ? Pose2D.OLIGO_MODE_DIMER : tc['fold_mode'];
    }

    public get_oligo_name(): string {
        let tc: any = this.get_target_conditions();
        if (tc == null) return null;
        return tc.hasOwnProperty('oligo_name') ? tc['oligo_name'] : null;
    }

    public set_target_oligos(target_oligos: any[]): void {
        this._target_oligos = target_oligos == null ? null : JSON.parse(JSON.stringify(target_oligos));
    }

    public set_target_oligo(target_oligo: any[]): void {
        this._target_oligo = target_oligo == null ? null : target_oligo.slice();
    }

    public get_oligo_order(): any[] {
        return this._oligo_order;
    }

    public set_oligo_order(oligo_order: any[]): void {
        this._oligo_order = oligo_order == null ? null : oligo_order.slice();
    }

    public get_oligos_paired(): number {
        return this._oligos_paired;
    }

    public set_oligos_paired(oligos_paired: number): void {
        this._oligos_paired = oligos_paired;
    }

    public get_target_pairs(): any[] {
        return this._target_pairs;
    }

    public set_target_pairs(target_pairs: any[]): void {
        this._target_pairs = target_pairs.slice();
    }

    public get_target_oligo_order(): any[] {
        return this._target_oligo_order;
    }

    public set_target_oligo_order(oligo_order: any[]): void {
        this._target_oligo_order = oligo_order == null ? null : oligo_order.slice();
    }

    public get_sequence(): number[] {
        return this._sequence;
    }

    public set_sequence(seq: number[]): void {
        this._sequence = seq.slice();
    }

    public get_puzzle_locks(): boolean[] {
        return this._puzzle_locks;
    }

    public set_puzzle_locks(locks: boolean[]): void {
        this._puzzle_locks = locks;
    }

    public get_forced_struct(): any[] {
        return this._forced_struct;
    }

    public set_forced_struct(forced: any[]): void {
        this._forced_struct = forced;
    }

    public set_target_conditions(target_conditions: any): void {
        this._target_conditions = JSON.stringify(target_conditions);
    }

    public get_target_conditions(): any {
        return (this._target_conditions == null ? null : JSON.parse(this._target_conditions));
    }

    public get_pairs(temp: number = 37): number[] {
        return this._pairs_array[temp];
    }

    public get_stable(): boolean {
        return this._stable;
    }

    public get_param(index: UndoBlockParam, temp: number = 37): any {
        if (this._params_array[temp] != null) {
            return this._params_array[temp][index];
        } else {
            return null;
        }
    }

    public set_pairs(pairs: number[], temp: number = 37): void {
        this._pairs_array[temp] = pairs.slice();
    }

    public set_stable(stable: boolean): void {
        this._stable = stable;
    }

    public set_param(index: UndoBlockParam, val: any, temp: number = 37): void {
        if (this._params_array[temp] == null) {
            this._params_array[temp] = [];
        }
        this._params_array[temp][index] = val;
    }

    public set_basics(folder: Folder, temp: number = 37): void {
        let best_pairs: any[] = this.get_pairs(temp);
        let seq: any[] = this._sequence;

        this.set_param(UndoBlockParam.PARAM_GU, EPars.num_gu_pairs(seq, best_pairs), temp);
        this.set_param(UndoBlockParam.PARAM_GC, EPars.num_gc_pairs(seq, best_pairs), temp);
        this.set_param(UndoBlockParam.PARAM_AU, EPars.num_ua_pairs(seq, best_pairs), temp);
        this.set_param(UndoBlockParam.PARAM_STACK, EPars.get_longest_stack_length(best_pairs), temp);
        this.set_param(UndoBlockParam.PARAM_REPETITION, EPars.get_sequence_repetition(EPars.sequence_array_to_string(seq), 5), temp);
        let full_seq: any[] = seq.slice();
        if (this._target_oligo) {
            if (this.get_oligo_mode() == Pose2D.OLIGO_MODE_DIMER) full_seq.push(EPars.RNABASE_CUT);
            if (this.get_oligo_mode() == Pose2D.OLIGO_MODE_EXT5P) {
                full_seq = this._target_oligo.concat(full_seq);
            } else {
                full_seq = full_seq.concat(this._target_oligo);
            }
        } else if (this._target_oligos) {
            for (let ii: number = 0; ii < this._target_oligos.length; ii++) {
                full_seq.push(EPars.RNABASE_CUT);
                full_seq = full_seq.concat(this._target_oligos[this._oligo_order[ii]].sequence);
            }
        }
        let nnfe: number[] = [];
        let total_fe: number = folder.score_structures(full_seq, best_pairs, temp, nnfe);
        this.set_param(UndoBlockParam.PARAM_FE, total_fe, temp);
        this.set_param(UndoBlockParam.PARAM_NNFE_ARRAY, nnfe, temp);
    }

    public set_meltingpoint_and_dotplot(folder: Folder): void {
        let pose_seq: string = EPars.sequence_array_to_string(this._sequence);

        let datablock: UndoBlock = this;

        if (datablock.get_param(UndoBlockParam.PARAM_DOTPLOT, 37) == null) {
            let dot_array: any[] = folder.get_dot_plot(datablock.get_sequence(), datablock.get_pairs(37), 37);
            datablock.set_param(UndoBlockParam.PARAM_DOTPLOT, dot_array, 37);
            this._dotplot.set_type(PlotType.SCATTER);
            this._dotplot.set_2d_data(dot_array, pose_seq.length);
        }

        for (let ii = 37; ii < 100; ii += 10) {
            if (datablock.get_pairs(ii) == null) {
                datablock.set_pairs(folder.fold_sequence(datablock.get_sequence(), null, null, ii), ii);
            }

            if (datablock.get_param(UndoBlockParam.PARAM_DOTPLOT, ii) == null) {
                let dot_temp_array: any[] = folder.get_dot_plot(datablock.get_sequence(), datablock.get_pairs(ii), ii);
                datablock.set_param(UndoBlockParam.PARAM_DOTPLOT, dot_temp_array, ii);
            }
        }

        let ref_pairs: number[] = datablock.get_pairs(37);

        let pair_scores: number[] = [];
        let max_pair_scores: number[] = [];

        for (let ii = 37; ii < 100; ii += 10) {
            if (datablock.get_param(UndoBlockParam.PARAM_PROB_SCORE, ii)) {
                pair_scores.push(1 - datablock.get_param(UndoBlockParam.PARAM_PAIR_SCORE, ii));
                max_pair_scores.push(1.0);
                continue;
            }
            let cur_dat: number[] = datablock.get_param(UndoBlockParam.PARAM_DOTPLOT, ii);
            let cur_pairs: number[] = datablock.get_pairs(ii);
            let prob_score: number = 0;
            let score_count: number = 0;

            for (let jj: number = 0; jj < cur_dat.length; jj += 3) {
                let index_i: number = cur_dat[jj] - 1;
                let index_j: number = cur_dat[jj + 1] - 1;

                if (index_i < index_j) {
                    if (ref_pairs[index_i] == index_j) {
                        prob_score += Number(cur_dat[jj + 2]);
                        score_count++;
                    }
                } else if (index_j < index_i) {
                    if (ref_pairs[index_j] == index_i) {
                        prob_score += Number(cur_dat[jj + 2]);
                        score_count++;
                    }
                }
            }

            if (score_count > 0) {
                prob_score /= score_count;
            }

            let num_paired: number = 0;
            for (let jj = 0; jj < cur_pairs.length; jj++) {
                if (cur_pairs[jj] > jj) {
                    num_paired += 2;
                }
            }
            let pair_score: number = Number(num_paired) / ref_pairs.length;

            pair_scores.push(1 - pair_score);
            max_pair_scores.push(1.0);

            datablock.set_param(UndoBlockParam.PARAM_PROB_SCORE, prob_score, ii);
            datablock.set_param(UndoBlockParam.PARAM_PAIR_SCORE, pair_score, ii);
        }

        this._meltplot.set_type(PlotType.LINE);
        this._meltplot.set_data(pair_scores, max_pair_scores);

        let init_score: number = datablock.get_param(UndoBlockParam.PARAM_PROB_SCORE, 37);

        let meltpoint: number = 107;
        for (let ii = 47; ii < 100; ii += 10) {
            let current_score: number = datablock.get_param(UndoBlockParam.PARAM_PROB_SCORE, ii);
            if (current_score < init_score * 0.5) {
                meltpoint = ii;
                break;
            }
        }

        datablock.set_param(UndoBlockParam.PARAM_MELTING_POINT, meltpoint, 37);
    }

    public get_dotplot(): Plot {
        return this._dotplot;
    }

    public get_meltplot(): Plot {
        return this._meltplot;
    }

    public get_order_map(other_order: any[]): any[] {
        if (this._target_oligos == null) return null;

        let idx_map: any[] = [];
        let ofs: any[] = [];
        let ii: number = this._sequence.length;
        for (let jj = 0; jj < this._target_oligos.length; jj++) {
            let kk = (other_order == null ? jj : other_order[jj]);
            ofs[kk] = ii;
            ii += 1 + this._target_oligos[kk].sequence.length;
        }
        for (let ii = 0; ii < this._sequence.length; ii++) idx_map[ii] = ii;
        for (let jj = 0; jj < this._target_oligos.length; jj++) {
            let kk = ofs[jj];
            let xx: number;
            for (xx = 0; xx <= this._target_oligos[jj].sequence.length; xx++) {
                idx_map[ii + xx] = kk + xx;
            }
            ii += xx;
        }
        return idx_map;
    }

    private _sequence: number[];
    private _pairs_array: number[][] = [];
    private _params_array: number[][] = [];
    private _stable: boolean = false;
    private _target_oligo: any[] = null;
    private _target_oligos: any[] = null;
    private _oligo_order: any[] = null;
    private _oligos_paired: number = 0;
    private _target_pairs: number[] = [];
    private _target_oligo_order: any[] = null;
    private _puzzle_locks: boolean[] = [];
    private _forced_struct: any[] = [];
    private _target_conditions: string = null;
    private _dotplot: Plot = new Plot();
    private _meltplot: Plot = new Plot();
}