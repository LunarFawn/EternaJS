import {Point, Sprite, Texture} from "pixi.js";
import {Updatable} from "../../flashbang/core/Updatable";
import {Vector2} from "../../flashbang/geom/Vector2";
import {ContainerObject} from "../../flashbang/objects/ContainerObject";
import {Pose2D} from "../pose2D/Pose2D";
import {BitmapManager} from "../resources/BitmapManager";
import {Bitmaps} from "../resources/Bitmaps";

export class NovaPaintHint extends ContainerObject implements Updatable {
    constructor(start: Point, end: Point, loop: boolean) {
        super();
        this._active = false;
        this._startPoint = start;
        this._endPoint = end;
        this._loop = loop;

        this._no_click = BitmapManager.get_bitmap(Bitmaps.NovaFinger);
        this._click_img = BitmapManager.get_bitmap(Bitmaps.NovaFingerClick);

        this._img = new Sprite(this._no_click);
        this.container.addChild(this._img);
    }

    public InitializeHint(): void {
        this._active = true;
        this._totalDistance = 0.0;
        this.display.position = this._startPoint;
    }

    public set_anchor_nucleotide(rna: Pose2D, base: number): void {
        this._rna = rna;
        this._base = base;
        this._anchor_set = true;
    }

    public update(dt :number): void {
        if (!this._active) {
            return;
        }

        let current_time = this.mode.time;

        let startPos: Point = this._startPoint;
        if (this._anchor_set) {
            startPos = this._rna.get_base_xy(this._base);
        }

        if (this._lastTimeTick == 0) {
            this._lastTimeTick = current_time;
            return;
        }

        if (this._startAnimTime == -1) {
            this._startAnimTime = current_time;
        }

        let stageTime: number = (current_time - this._startAnimTime);
        let dir: Vector2 = new Vector2(this._endPoint.x - this._startPoint.x, this._endPoint.y - this._startPoint.y);
        if (stageTime < 1.5 && this._curStage == 0) {
            if (stageTime >= 1.4) {
                ++this._curStage;
            } else if (stageTime > 0.7) {
                this._img.texture = this._click_img;
            }
        } else if (this._curStage == 1) {
            let deltaTime: number = current_time - this._lastTimeTick;
            // Move from our current position to the end
            let stepDistance = deltaTime * NovaPaintHint.PAINT_HINT_SPEED;
            this._totalDistance += stepDistance;
            if (this._totalDistance >= Vector2.distance(this._startPoint.x, this._startPoint.y, this._endPoint.x, this._endPoint.y) - 1.5) {
                this._endAnimTime = current_time;
                ++this._curStage;
            }
        } else if (this._curStage == 2) {
            if (!this._loop) {
                this._active = false;
                return;
            }

            let endTime: number = (current_time - this._endAnimTime) / 1000.0;
            if (endTime > 1.0) {
                this._startAnimTime = -1;
                this._curStage = 0;
                this.InitializeHint();
            } else if (endTime > 0.5) {
                this._img.texture = this._no_click;
            }
        }
        this._lastTimeTick = current_time;

        dir.length = this._totalDistance;
        this.display.position = new Point(startPos.x + dir.x, startPos.y + dir.y);
    }

    private readonly _startPoint: Point;
    private readonly _loop: boolean;
    private readonly _img: Sprite;
    private readonly _no_click: Texture;
    private readonly _click_img: Texture;

    private _rna: Pose2D;
    private _base: number;
    private _anchor_set: boolean = false;
    private _endPoint: Point;
    private _active: boolean;
    private _lastTimeTick: number = 0;
    private _startAnimTime: number = -1;
    private _endAnimTime: number = -1;
    private _curStage: number = 0;
    private _totalDistance: number = 0;

    private static readonly PAINT_HINT_SPEED: number = 80;
}
