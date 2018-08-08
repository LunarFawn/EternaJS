import {ROPWait} from "./ROPWait";
import {RScriptOp} from "./RScriptOp";

/*
 * RScript Operation Tree.
 * Clearly organizes all of the input operations into a tree.
 * :: Note - Use a tree so that we can support if conditionals easily in the near future.
 * :: Note - Only one entry point.
 */
export class RScriptOpTree {
    public AddNode(node: RScriptOp): void {
        if (!node) {
            return;
        }

        if (this._head == null) {
            this._head = node;
            this._curptr = node;
        } else {
            this._curptr.AddChildOp(node);
            this._curptr = node;
        }
    }

    public FinishCreation(): void {
        this._curptr = this._head;
    }

    public next(): RScriptOp {
        if (!this._curptr) {
            return null;
        }

        if (this._curptr instanceof ROPWait) {
            this._curptr.exec();
            let waitRet: RScriptOp = this._curptr.get_pause_next();
            if (waitRet !== this._curptr && this._curptr.IsPaused() && waitRet instanceof ROPWait) {
                // If the next instruction can be executed (as determined by ROPWait),
                // then execute it.
                this._waitQueue.push(this._curptr);
                this._curptr = waitRet;
                return waitRet;
            } else if (this._curptr.IsPaused() && this._waitQueue.indexOf(this._curptr)) {
                this._waitQueue.push(this._curptr);
                return null;
            } else {
                // If it cannot then see if the wait queue is clear.
                if (this.CheckWaitQueueContinue()) {
                    // Clear queue and proceed.
                    this._waitQueue.splice(0);
                } else {
                    return null;
                }
            }
        }

        let ret: RScriptOp = this._curptr;
        this._curptr = this._curptr.next();
        return ret;
    }

    private CheckWaitQueueContinue(): boolean {
        for (let op of this._waitQueue) {
            if (op.IsPaused()) {
                return false;
            }
        }
        return true;
    }

    private _head: RScriptOp;
    private _curptr: RScriptOp;
    private _waitQueue: RScriptOp[] = [];
}