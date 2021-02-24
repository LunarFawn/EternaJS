import * as log from 'loglevel';
import {
    Container, Graphics, Point, Sprite, Texture, Rectangle
} from 'pixi.js';
import {Registration, Value} from 'signals';
import EPars, {RNABase, RNAPaint} from 'eterna/EPars';
import Eterna from 'eterna/Eterna';
import ExpPainter from 'eterna/ExpPainter';
import {
    ContainerObject, InputUtil, Flashbang, Dragger, DisplayUtil, SceneObject, SerialTask, Easing,
    ParallelTask, AlphaTask, LocationTask, DelayTask, SelfDestructTask, Vector2, Arrays,
    RepeatingTask, Updatable, Assert
} from 'flashbang';
import {Move} from 'eterna/mode/PoseEdit/PoseEditMode';
import LightRay from 'eterna/vfx/LightRay';
import TextBalloon from 'eterna/ui/TextBalloon';
import ROPWait from 'eterna/rscript/ROPWait';
import Fonts from 'eterna/util/Fonts';
import Sounds from 'eterna/resources/Sounds';
import BaseGlow from 'eterna/vfx/BaseGlow';
import BitmapManager from 'eterna/resources/BitmapManager';
import Booster from 'eterna/mode/PoseEdit/Booster';
import GameMode from 'eterna/mode/GameMode';
import Utility from 'eterna/util/Utility';
import Folder from 'eterna/folding/Folder';
import SecStruct from 'eterna/rnatypes/SecStruct';
import Sequence from 'eterna/rnatypes/Sequence';
import Base from './Base';
import BaseDrawFlags from './BaseDrawFlags';
import EnergyScoreDisplay from './EnergyScoreDisplay';
import HighlightBox, {HighlightType} from './HighlightBox';
import BaseRope from './BaseRope';
import PseudoknotLines from './PseudoknotLines';
import Molecule from './Molecule';
import PaintCursor from './PaintCursor';
import PoseField from './PoseField';
import PoseUtil from './PoseUtil';
import PuzzleEditOp from './PuzzleEditOp';
import RNAAnchorObject from './RNAAnchorObject';
import RNALayout, {RNATreeNode} from './RNALayout';
import ScoreDisplayNode, {ScoreDisplayNodeType} from './ScoreDisplayNode';
import ExplosionFactorPanel from './ExplosionFactorPanel';
import triangulate from './triangulate';
import Tooltips from '../ui/Tooltips';
import AnnotationItem, {
    AnnotationData,
    AnnotationDataCollection,
    AnnotationCategory,
    AnnotationRange,
    AnnotationItemType,
    AnnotationDisplayObject,
    AnnotationGraph,
    AnnotationGraphNode,
    AnnotationPlacement,
    AnnotationPosition,
    AnnotationBaseConflict,
    AnnotationPositionConflict,
    AnnotationArguments,
    AnnotationBaseConflicts
} from '../ui/AnnotationItem';
import AnnotationCard from '../ui/AnnotationCard';

type InteractionEvent = PIXI.InteractionEvent;

interface Mut {
    pos: number;
    base: string;
}

export enum Layout {
    MOVE,
    ROTATE_STEM,
    FLIP_STEM
}

export type PoseMouseDownCallback = (e: InteractionEvent, closestDist: number, closestIndex: number) => void;

export default class Pose2D extends ContainerObject implements Updatable {
    public static readonly COLOR_CURSOR: number = 0xFFC0CB;
    public static readonly ZOOM_SPACINGS: number[] = [45, 30, 20, 14, 7];

    public static readonly OLIGO_MODE_DIMER: number = 1;
    public static readonly OLIGO_MODE_EXT3P: number = 2;
    public static readonly OLIGO_MODE_EXT5P: number = 3;

    private static readonly SCORES_POSITION_Y = 128;

    public readonly onCreateAnnotation: Value<AnnotationArguments> = new Value<AnnotationArguments>({ranges: []});
    public readonly onSelectAnnotation: Value<AnnotationData> = new Value<AnnotationData>({
        id: '',
        type: AnnotationItemType.ANNOTATION,
        timestamp: 0,
        playerID: Eterna.playerID,
        title: '',
        ranges: []
    });

    public readonly onSelectLayer: Value<AnnotationData> = new Value<AnnotationData>({
        id: '',
        type: AnnotationItemType.LAYER,
        title: '',
        playerID: Eterna.playerID
    });

    public readonly onEditAnnotation: Value<AnnotationData | null> = new Value<AnnotationData | null>(null);

    constructor(poseField: PoseField, editable: boolean) {
        super();
        this._poseField = poseField;
        this._editable = editable;
    }

    protected added() {
        super.added();

        this._scoreNodeHighlight = new Graphics();
        this.container.addChild(this._scoreNodeHighlight);

        this._baseRope = new BaseRope(this);
        this.addObject(this._baseRope, this.container);

        this._pseudoknotLines = new PseudoknotLines(this);
        this.addObject(this._pseudoknotLines, this.container);

        this.container.addChild(this._baseLayer);

        this._primaryScoreEnergyDisplay = new EnergyScoreDisplay(111, 40);
        this._primaryScoreEnergyDisplay.position = new Point(17, Pose2D.SCORES_POSITION_Y);
        this.container.addChild(this._primaryScoreEnergyDisplay);

        this._deltaScoreEnergyDisplay = new EnergyScoreDisplay(111, 40);
        this._deltaScoreEnergyDisplay.position = new Point(17 + 119, Pose2D.SCORES_POSITION_Y);
        this._deltaScoreEnergyDisplay.visible = false;
        this.container.addChild(this._deltaScoreEnergyDisplay);

        this._secondaryScoreEnergyDisplay = new EnergyScoreDisplay(111, 40);
        this._secondaryScoreEnergyDisplay.position = new Point(17 + 119 * 2, Pose2D.SCORES_POSITION_Y);
        this._secondaryScoreEnergyDisplay.visible = false;
        this.container.addChild(this._secondaryScoreEnergyDisplay);

        this._moleculeLayer = new Container();
        this.container.addChild(this._moleculeLayer);
        this._moleculeLayer.visible = false;

        this._energyTextLayer = new Container();
        this.container.addChild(this._energyTextLayer);

        this._paintCursor = new PaintCursor();
        this._paintCursor.display.visible = false;
        this.addObject(this._paintCursor, this.container);

        this._explosionFactorPanel = new ExplosionFactorPanel();
        this._explosionFactorPanel.display.position = new Point(17, Pose2D.SCORES_POSITION_Y + 82);
        this._explosionFactorPanel.display.visible = false;
        this._explosionFactorPanel.factorUpdated.connect((factor: number) => {
            this._explosionFactor = factor;
            this.computeLayout(true);
            this._redraw = true;
        });
        this.addObject(this._explosionFactorPanel, this.container);

        this._explosionRays = [];
        for (let ii = 0; ii < 10; ii++) {
            const ray = new LightRay();
            this._explosionRays.push(ray);
            this.addObject(ray, this.container);
        }

        this._selectionHighlightBox = new HighlightBox(this, HighlightType.DESIGN);
        this.addObject(this._selectionHighlightBox, this.container);

        this._restrictedHighlightBox = new HighlightBox(this, HighlightType.RESTRICTED);
        this.addObject(this._restrictedHighlightBox, this.container);

        this._unstableHighlightBox = new HighlightBox(this, HighlightType.UNSTABLE);
        this.addObject(this._unstableHighlightBox, this.container);

        this._userDefinedHighlightBox = new HighlightBox(this, HighlightType.USER_DEFINED);
        this.addObject(this._userDefinedHighlightBox, this.container);

        this._forcedHighlightBox = new HighlightBox(this, HighlightType.FORCED);
        this.addObject(this._forcedHighlightBox, this.container);

        this._shiftHighlightBox = new HighlightBox(this, HighlightType.SHIFT);
        this.addObject(this._shiftHighlightBox, this.container);

        this._annotationHighlightBox = new HighlightBox(this, HighlightType.ANNOTATION);
        this._annotationHighlightBox.display.cursor = 'pointer';
        this.addObject(this._annotationHighlightBox, this.container);
        if (Tooltips.instance != null) {
            this.regs.add(Tooltips.instance.addTooltip(this._annotationHighlightBox, 'Create Annotation'));
        }

        if (!this._editable) {
            this._currentColor = -1;
        }

        this.annotationCanvas = new Graphics();
        this.container.addChild(this.annotationCanvas);

        this._strandLabel = new TextBalloon('', 0x0, 0.8);
        this._strandLabel.display.visible = false;
        this.addObject(this._strandLabel, this.container);

        this.pointerMove.connect((p) => this.onMouseMoved(p.data.global));
        this.pointerDown.filter(InputUtil.IsLeftMouse).connect((e) => this.callStartMousedownCallback(e));
        this.pointerOut.connect(() => this.onMouseOut());

        // handle view settings
        this.regs.add(Eterna.settings.showNumbers.connectNotify((value) => { this.showNumbering = value; }));
        this.regs.add(Eterna.settings.showRope.connectNotify((value) => { this.showRope = value; }));
        this.regs.add(Eterna.settings.showLetters.connectNotify((value) => { this.lettermode = value; }));
        this.regs.add(
            Eterna.settings.useContinuousColors.connectNotify((value) => { this.useContinuousExpColors = value; })
        );
        this.regs.add(Eterna.settings.useExtendedColors.connectNotify((value) => { this.useExtendedScale = value; }));
        this.regs.add(
            Eterna.settings.displayFreeEnergies.connectNotify((value) => { this.displayScoreTexts = value; })
        );
        this.regs.add(
            Eterna.settings.highlightRestricted.connectNotify((value) => { this.highlightRestricted = value; })
        );
        this.regs.add(Eterna.settings.simpleGraphics.connectNotify((value) => { this.useSimpleGraphics = value; }));
        this.regs.add(Eterna.settings.annotationModeActive.connectNotify((value) => {
            this.annotationModeActive = value;
        }));
        this.regs.add(Eterna.settings.usePuzzlerLayout.connect(() => this.computeLayout()));
    }

    public setSize(width: number, height: number): void {
        this._width = width;
        this._height = height;

        this.container.hitArea = new Rectangle(0, 0, width, height);
    }

    public get primaryScoreDisplay(): EnergyScoreDisplay {
        return this._primaryScoreEnergyDisplay;
    }

    public get secondaryScoreDisplay(): EnergyScoreDisplay {
        return this._secondaryScoreEnergyDisplay;
    }

    public addAnchoredObject(obj: RNAAnchorObject): void {
        this._anchoredObjects.push(obj);
    }

    public removeAnchoredObject(obj: RNAAnchorObject): void {
        for (let ii = 0; ii < this._anchoredObjects.length; ++ii) {
            if (obj === this._anchoredObjects[ii]) {
                this._anchoredObjects.splice(ii, 1);
                break;
            }
        }
    }

    public get isAnimating(): boolean {
        return this._baseToX != null;
    }

    public get isFolding(): boolean {
        return (this.lastSampledTime - this._foldStartTime < this._foldDuration);
    }

    public get annotations(): AnnotationDisplayObject[] {
        return this._annotations;
    }

    public get layers(): AnnotationDisplayObject[] {
        return this._layers;
    }

    public visualizeFeedback(dat: number[], mid: number, lo: number, hi: number, startIndex: number): void {
        // coloring
        const newdat: number[] = ExpPainter.transformData(dat, hi, lo);
        this._expPainter = new ExpPainter(newdat, startIndex);
        this._expMid = mid;
        this._expHi = hi;
        this.paintFeedback();
    }

    public paintFeedback(): void {
        if (!this._expPainter) {
            return;
        }

        this._expPainter.continuous = this._expContinuous;
        this._expPainter.extendedScale = this._expExtendedScale;

        for (let ii = 0; ii < this._sequence.length; ii++) {
            this._bases[ii].setColorLevel(
                true, this._expPainter.getColorLevelWithMidpoint(ii, this._expMid, this._expHi)
            );
        }
        this._redraw = true;
    }

    public clearFeedback(): void {
        for (let ii = 0; ii < this._sequence.length; ii++) {
            this._bases[ii].setColorLevel(false, -1);
        }
        this._redraw = true;
    }

    public get zoomLevel(): number {
        return this._zoomLevel;
    }

    public setZoomLevel(zoomLevel: number, animate: boolean = true, center: boolean = false): void {
        if ((this._zoomLevel !== zoomLevel || center) && animate) {
            if (this._zoomLevel === zoomLevel && center) {
                if (Math.abs(this._width / 2 - this._offX) + Math.abs(this._height / 2 - this._offY) < 50) {
                    return;
                }
            }

            // Update Annotations
            if (this._annotations.length > 0) {
                this.updateAnnotationSpaceAvailability();
                this.eraseAnnotations(true);
                this.drawAnnotations();
            }

            this._startOffsetX = this._offX;
            this._startOffsetY = this._offY;

            let scaler = 1;
            if (zoomLevel > this._zoomLevel) {
                scaler = Pose2D.ZOOM_SPACINGS[zoomLevel] / Pose2D.ZOOM_SPACINGS[this._zoomLevel];
            }

            if (!this._offsetTranslating && !center) {
                this._endOffsetX = scaler * (this._offX - this._width / 2) + this._width / 2;
                this._endOffsetY = scaler * (this._offY - this._height / 2) + this._height / 2;
            } else if (this._offsetTranslating) {
                this._endOffsetX = scaler * (this._endOffsetX - this._width / 2) + this._width / 2;
                this._endOffsetY = scaler * (this._endOffsetY - this._height / 2) + this._height / 2;
            } else {
                this._endOffsetX = this._width / 2;
                this._endOffsetY = this._height / 2;
            }

            this._offsetTranslating = true;

            this._zoomLevel = zoomLevel;
            this.computeLayout(true);
            this._redraw = true;
        } else if (this._zoomLevel !== zoomLevel) {
            this._zoomLevel = zoomLevel;
            this.computeLayout(true);
            this._redraw = true;
        }
    }

    public computeDefaultZoomLevel(): number {
        const n: number = this.fullSequenceLength;
        const xarray: number[] = new Array(n);
        const yarray: number[] = new Array(n);

        const rnaCoords: RNALayout = new RNALayout(Pose2D.ZOOM_SPACINGS[0], Pose2D.ZOOM_SPACINGS[0]);
        rnaCoords.setupTree(this._pairs, this._targetPairs);
        rnaCoords.drawTree(this._customLayout);
        rnaCoords.getCoords(xarray, yarray);

        const xmin: number = Math.min(...xarray);// xarray[0];
        const xmax: number = Math.max(...xarray);
        const ymin: number = Math.min(...yarray);
        const ymax: number = Math.max(...yarray);
        const xdiff: number = xmax - xmin;
        const ydiff: number = ymax - ymin;
        const xscale: number = xdiff / this._width;
        const yscale: number = ydiff / this._height;

        const scale: number = Math.max(xscale, yscale);
        if (scale < 1.0) {
            return 0;
        } else if ((30 / 45) * scale < 1.0) {
            return 1;
        } else if ((20 / 45) * scale < 1.0) {
            return 2;
        } else if ((14 / 45) * scale < 1.0) {
            return 3;
        } else {
            return 4;
        }
    }

    public set currentColor(col: RNAPaint) {
        this._currentColor = col;
    }

    public get currentColor(): RNAPaint {
        return this._currentColor;
    }

    public set currentArrangementTool(col: Layout) {
        this._currentArrangementTool = col;
    }

    public get currentArrangementTool(): Layout {
        return this._currentArrangementTool;
    }

    public doneColoring(): void {
        this._coloring = false;

        let needUpdate = false;

        if (this._customLayoutChanged) {
            this.checkPairs();
            this.updateMolecule();
            this.generateScoreNodes();
            this.callPoseEditCallback();

            // Update Annotations
            if (this._annotations.length > 0) {
                this.eraseAnnotations(false);
                this.drawAnnotations();
            }
            return;
        }

        if (this._mutatedSequence == null) {
            return;
        }

        if (this._mutatedSequence.length !== this.fullSequenceLength) {
            throw new Error("Mutated sequence and original sequence lengths don't match");
        }

        let div = 1;
        if (this._currentColor === RNAPaint.PAIR
            || this._currentColor === RNAPaint.GC_PAIR
            || this._currentColor === RNAPaint.AU_PAIR
            || this._currentColor === RNAPaint.GU_PAIR) {
            div = 2;
        }

        const offset: number = (
            this._oligo != null
            && this._oligoMode === Pose2D.OLIGO_MODE_EXT5P
        ) ? this._oligo.length : 0;
        let numMut = 0;
        const muts: Mut[] = [];
        for (let ii = 0; ii < this._sequence.length; ii++) {
            if (this._sequence.nt(ii) !== this._mutatedSequence.nt(ii + offset)) {
                numMut++;
                this._sequence.setNt(ii, this._mutatedSequence.nt(ii + offset));
                muts.push({pos: ii + 1, base: EPars.nucleotideToString(this._sequence.nt(ii))});
                needUpdate = true;
            }
        }
        if (needUpdate) {
            this.callTrackMovesCallback(numMut / div, muts);
        }
        if (
            needUpdate
            || this._lockUpdated
            || this._bindingSiteUpdated
            || this._designStructUpdated
        ) {
            this.checkPairs();
            this.updateMolecule();
            this.generateScoreNodes();
            this.callPoseEditCallback();

            // Update Annotations
            if (this._annotations.length > 0) {
                this.eraseAnnotations();
                this.drawAnnotations();
            }
        }

        this._mutatedSequence = null;
        this._lockUpdated = false;
        this._bindingSiteUpdated = false;
        this._designStructUpdated = false;
    }

    public setMutated(seqArr: Sequence): void {
        Assert.assertIsDefined(this._mutatedSequence);
        const n: number = Math.min(this._mutatedSequence.length, seqArr.length);
        const offset: number = (
            this._oligo != null && this._oligoMode === Pose2D.OLIGO_MODE_EXT5P
        ) ? this._oligo.length : 0;

        for (let ii = 0; ii < n; ii++) {
            if (this._mutatedSequence.nt(ii) !== seqArr.nt(ii) && !this.isLocked(offset + ii)) {
                this._mutatedSequence.setNt(ii, seqArr.nt(ii));
                this._bases[offset + ii].setType(seqArr.nt(ii));
            }
        }
    }

    public pasteSequence(sequence: Sequence): void {
        if (sequence == null) {
            return;
        }

        let numMut = 0;
        const muts: Mut[] = [];

        const n: number = Math.min(sequence.length, this._sequence.length);
        let needUpdate = false;
        const offset: number = (
            this._oligo != null && this._oligoMode === Pose2D.OLIGO_MODE_EXT5P
        ) ? this._oligo.length : 0;

        for (let ii = 0; ii < n; ii++) {
            if (sequence.nt(ii) === RNABase.UNDEFINED) continue;
            if (this._sequence.nt(ii) !== sequence.nt(ii) && !this.isLocked(offset + ii)) {
                numMut++;
                this._sequence.setNt(ii, sequence.nt(ii));
                muts.push({pos: ii + 1, base: EPars.nucleotideToString(this._sequence.nt(ii))});
                this._bases[offset + ii].setType(sequence.nt(ii));
                needUpdate = true;
            }
        }

        if (needUpdate) {
            this.callTrackMovesCallback(numMut, muts);

            this.checkPairs();
            this.updateMolecule();
            this.generateScoreNodes();
            this.callPoseEditCallback();

            // Update Annotations
            if (this._annotations.length > 0) {
                this.updateAnnotationSpaceAvailability();
                this.eraseAnnotations(true);
                this.drawAnnotations();
            }
        }
    }

    public getBaseLoc(seq: number, out: Point | null = null): Point {
        if (out == null) {
            out = new Point();
        }
        out.x = this._bases[seq].x + this._offX;
        out.y = this._bases[seq].y + this._offY;
        return out;
    }

    public getEnergyScorePos(index: number, out: Point | null = null): Point {
        if (out === null) {
            out = new Point();
        }
        Assert.assertIsDefined(
            this._scoreTexts,
            "Can't get substructure score position, because the scores do not exist"
        );
        out.x = this._scoreTexts[index].x;
        out.y = this._scoreTexts[index].y;
        return out;
    }

    public getBaseOutXY(seq: number, out: Point | null = null): Point {
        out = this._bases[seq].getOutXY(out);
        out.x += this._offX;
        out.y += this._offY;
        return out;
    }

    public clearMouse(): void {
        // document.getElementById(Eterna.PIXI_CONTAINER_ID).style.cursor = '';
        this._paintCursor.display.visible = false;
        this._strandLabel.display.visible = false;
    }

    public parseCommand(command: RNAPaint, closestIndex: number): [string, PuzzleEditOp, RNABase[]?] | null {
        switch (command) {
            case RNAPaint.ADD_BASE:
                return PoseUtil.addBaseWithIndex(closestIndex, this._pairs);

            case RNAPaint.ADD_PAIR:
                return PoseUtil.addPairWithIndex(closestIndex, this._pairs);

            case RNAPaint.DELETE:
                return this.deleteBaseWithIndex(closestIndex);

            default:
                return null;
        }
    }

    public parseCommandWithPairs(
        command: RNAPaint, closestIndex: number, pairs: SecStruct
    ): [string, PuzzleEditOp, RNABase[]?] | null {
        switch (command) {
            case RNAPaint.ADD_BASE:
                return PoseUtil.addBaseWithIndex(closestIndex, pairs);

            case RNAPaint.DELETE:
                return this.deleteBaseWithIndexPairs(closestIndex, pairs);

            default:
                return null;
        }
    }

    public onPoseMouseDownPropagate(e: InteractionEvent, closestIndex: number): void {
        const altDown: boolean = Flashbang.app.isAltKeyDown;
        const ctrlDown: boolean = Flashbang.app.isControlKeyDown || Flashbang.app.isMetaKeyDown;
        const ctrlDownOrBaseMarking = ctrlDown || this.currentColor === RNAPaint.BASE_MARK;

        if ((this._coloring && !altDown) || ctrlDownOrBaseMarking) {
            if (ctrlDownOrBaseMarking && closestIndex >= this.sequence.length) {
                return;
            }
            this.onPoseMouseDown(e, closestIndex);
        }
    }

    /**
     * Rotate the stem containing nucleotide idx. Save your results in the
     * customLayout. To achieve this: figure out what stem you're in (helper);
     * find its center (for an axis of rotation); figure out what orientation
     * is clockwise of its current orientation, then get there. We do not store
     * persistent stem orientations, because we would then have to reset them
     * elsewhere with every refold.
     *
     * @param idx
     */
    public rotateStem(startIdx: number): void {
        // If this idx is not paired, it won't be in a stem; return.
        if (!this._targetPairs.isPaired(startIdx)) {
            return;
        }

        // 1. Get coords and set up a customLayout
        const rnaCoords: RNALayout = new RNALayout(
            Pose2D.ZOOM_SPACINGS[this._zoomLevel], Pose2D.ZOOM_SPACINGS[this._zoomLevel]
        );
        // rnaCoords.setupTree(this._pairs.filterForPseudoknots(), this._targetPairs.filterForPseudoknots());
        rnaCoords.setupTree(this._pairs, this._targetPairs);
        rnaCoords.drawTree(this._customLayout);
        const xarray: number[] = new Array(this._bases.length);
        const yarray: number[] = new Array(this._bases.length);
        rnaCoords.getCoords(xarray, yarray);

        const localCustomLayout: ([number, number] | [null, null])[] = [];
        for (let ii = 0; ii < this._bases.length; ++ii) {
            if (xarray[ii] === undefined || yarray[ii] === undefined) continue;
            localCustomLayout.push([
                xarray[ii],
                yarray[ii]
            ]);
        }

        // id stem
        const stem = this._targetPairs.stemWith(startIdx);

        // What is center of stem? Average coordinate. Could simplify calculation
        // a little by finding the bases of median index first or something, but
        // that itself takes a bit of work. Unlikely to become bottleneck.
        const center = ((s: [number, number][]) => {
            let x = 0;
            let y = 0;
            for (const bp of s) {
                for (const idx of bp) {
                    x += this._bases[idx].x;
                    y += this._bases[idx].y;
                }
            }
            return [x / (s.length * 2), y / (s.length * 2)];
        })(stem);

        // Determine stem orientation. Really we only care about the "next"
        // orientation, the one we want to impose. We do this by orienting the
        // single bp that is kind of a minimal stem, guaranteed to be there.
        // If the vector from smaller to larger index has +x and -y, we rotate
        // to be +x0y. That goes to 0x+y, to -x0y, to 0x-y.
        const vecBP = ((s: [number, number][], n: number) => {
            for (const bp of s) {
                if (bp[0] === n || bp[1] === n) {
                    if (bp[0] < bp[1]) {
                        return [
                            this._bases[bp[1]].x - this._bases[bp[0]].x,
                            this._bases[bp[1]].y - this._bases[bp[0]].y
                        ];
                    } else {
                        return [
                            this._bases[bp[1]].x - this._bases[bp[0]].x,
                            this._bases[bp[1]].y - this._bases[bp[0]].y
                        ];
                    }
                }
            }
            return [0, 0];
        })(stem, startIdx);

        // Find and sort the (smallest, largest) bp. Let's assume it's first.
        const firstbp = stem[0][0] < stem[0][1]
            ? stem[0]
            : [stem[0][1], stem[0][0]];
        const pairSpace = Pose2D.ZOOM_SPACINGS[this._zoomLevel];
        const primSpace = Pose2D.ZOOM_SPACINGS[this._zoomLevel];

        // Calculate new stem positions
        if (vecBP[1] < 0) {
            // new orientation is bottom-to-top
            for (const bp of stem) {
                // First work with smaller value, which is either smaller or
                // bigger than the center bp
                this._bases[bp[0]].setXY(
                    center[0] - pairSpace / 2,
                    center[1] + (bp[0] - firstbp[0] + 0.5 - stem.length / 2) * primSpace
                );
                this._bases[bp[0]].setDirty();
                this._bases[bp[1]].setXY(
                    center[0] + pairSpace / 2,
                    center[1] + (bp[0] - firstbp[0] + 0.5 - stem.length / 2) * primSpace
                );
                this._bases[bp[1]].setDirty();
            }
        } else if (vecBP[0] > 0) {
            // new orientation is right-to-left
            for (const bp of stem) {
                // First work with smaller value, which is either smaller or
                // bigger than the center bp
                this._bases[bp[0]].setXY(
                    center[0] - (bp[0] - firstbp[0] + 0.5 - stem.length / 2) * primSpace,
                    center[1] - pairSpace / 2
                );
                this._bases[bp[0]].setDirty();
                this._bases[bp[1]].setXY(
                    center[0] - (bp[0] - firstbp[0] + 0.5 - stem.length / 2) * primSpace,
                    center[1] + pairSpace / 2
                );
                this._bases[bp[1]].setDirty();
            }
        } else if (vecBP[1] > 0) {
            // new orientation is top-to-bottom
            for (const bp of stem) {
                // First work with smaller value, which is either smaller or
                // bigger than the center bp
                this._bases[bp[0]].setXY(
                    center[0] + pairSpace / 2,
                    center[1] - (bp[0] - firstbp[0] + 0.5 - stem.length / 2) * primSpace
                );
                this._bases[bp[0]].setDirty();
                this._bases[bp[1]].setXY(
                    center[0] - pairSpace / 2,
                    center[1] - (bp[0] - firstbp[0] + 0.5 - stem.length / 2) * primSpace
                );
                this._bases[bp[1]].setDirty();
            }
        } else if (vecBP[0] < 0) {
            // new orientation is left-to-right
            for (const bp of stem) {
                // First work with smaller value, which is either smaller or
                // bigger than the center bp
                this._bases[bp[0]].setXY(
                    center[0] + (bp[0] - firstbp[0] + 0.5 - stem.length / 2) * primSpace,
                    center[1] + pairSpace / 2
                );
                this._bases[bp[0]].setDirty();
                this._bases[bp[1]].setXY(
                    center[0] + (bp[0] - firstbp[0] + 0.5 - stem.length / 2) * primSpace,
                    center[1] - pairSpace / 2
                );
                this._bases[bp[1]].setDirty();
            }
        }
        for (const bp of stem) {
            for (let ii = 0; ii < localCustomLayout.length; ++ii) {
                localCustomLayout[bp[0]] = [
                    localCustomLayout[ii][0] as number + this._bases[bp[0]].x - this._bases[ii].x,
                    localCustomLayout[ii][1] as number + this._bases[bp[0]].y - this._bases[ii].y
                ];
                localCustomLayout[bp[1]] = [
                    localCustomLayout[ii][0] as number + this._bases[bp[1]].x - this._bases[ii].x,
                    localCustomLayout[ii][1] as number + this._bases[bp[1]].y - this._bases[ii].y
                ];
            }
        }

        // Use setter to force redraw
        this.customLayout = localCustomLayout;
        this._baseRope.enabled = true;
        this._baseRope.redraw(true);
        this._pseudoknotLines.redraw(true);
        this._customLayoutChanged = true;
    }

    /**
     * Flip the stem containing nucleotide idx. Save your results in the
     * customLayout. To achieve this: swap the position of each nt with its bp
     * partner. We do not store persistent stem orientations, because we would
     * then have to reset them elsewhere with every refold.
     *
     * @param idx
     */
    public flipStem(startIdx: number): void {
        // If this idx is not paired, it won't be in a stem; return.
        if (!this._targetPairs.isPaired(startIdx)) {
            return;
        }

        // 1. Get coords and set up a customLayout
        const rnaCoords: RNALayout = new RNALayout(
            Pose2D.ZOOM_SPACINGS[this._zoomLevel], Pose2D.ZOOM_SPACINGS[this._zoomLevel]
        );
        rnaCoords.setupTree(this._pairs, this._targetPairs);
        rnaCoords.drawTree(this._customLayout);
        const xarray: number[] = new Array(this._bases.length);
        const yarray: number[] = new Array(this._bases.length);
        rnaCoords.getCoords(xarray, yarray);

        const localCustomLayout: ([number, number] | [null, null])[] = [];
        for (let ii = 0; ii < this._bases.length; ++ii) {
            if (xarray[ii] === undefined || yarray[ii] === undefined) continue;
            localCustomLayout.push([
                (xarray[ii]), // * (Pose2D.ZOOM_SPACINGS[0] / Pose2D.ZOOM_SPACINGS[this._zoomLevel]),
                (yarray[ii])// * (Pose2D.ZOOM_SPACINGS[0] / Pose2D.ZOOM_SPACINGS[this._zoomLevel])
            ]);
        }

        // id stem
        const stem = this._targetPairs.stemWith(startIdx);

        // Calculate new stem positions
        for (const bp of stem) {
            // First work with smaller value, which is either smaller or
            // bigger than the center bp
            const tmp = [
                this._bases[bp[0]].x,
                this._bases[bp[0]].y
            ];
            this._bases[bp[0]].setXY(
                this._bases[bp[1]].x,
                this._bases[bp[1]].y
            );
            this._bases[bp[0]].setDirty();
            this._bases[bp[1]].setXY(
                tmp[0],
                tmp[1]
            );
            this._bases[bp[1]].setDirty();
        }

        for (const bp of stem) {
            for (let ii = 0; ii < localCustomLayout.length; ++ii) {
                localCustomLayout[bp[0]] = [
                    localCustomLayout[ii][0] as number + this._bases[bp[0]].x - this._bases[ii].x,
                    localCustomLayout[ii][1] as number + this._bases[bp[0]].y - this._bases[ii].y
                ];
                localCustomLayout[bp[1]] = [
                    localCustomLayout[ii][0] as number + this._bases[bp[1]].x - this._bases[ii].x,
                    localCustomLayout[ii][1] as number + this._bases[bp[1]].y - this._bases[ii].y
                ];
            }
        }
        // Use setter to ensure structure update.
        this.customLayout = localCustomLayout;
        this._baseRope.enabled = true;
        this._baseRope.redraw(true);
        this._pseudoknotLines.redraw(true);
        this._customLayoutChanged = true;
    }

    /**
     * Snap every base to a grid of size pairSpace/5. This is a nice even number
     * such that you get SOME gradations that are smaller than the space between
     * paired bases, but not so much that the whole thing feels sloppy. Also
     * permits nice tetraloop layouts. Only affects paired nucleotides.
     *
     * @param idx
     */
    public snapToGrid(): void {
        const pairSpace = Pose2D.ZOOM_SPACINGS[this._zoomLevel];
        const gridSpace = pairSpace;

        // 1. Get coords and set up a customLayout
        const rnaCoords: RNALayout = new RNALayout(
            pairSpace, pairSpace
        );
        rnaCoords.setupTree(this._pairs, this._targetPairs);
        rnaCoords.drawTree(this._customLayout);
        const xarray: number[] = new Array(this._bases.length);
        const yarray: number[] = new Array(this._bases.length);
        rnaCoords.getCoords(xarray, yarray);

        const localCustomLayout: ([number, number] | [null, null])[] = [];
        for (let ii = 0; ii < this._bases.length; ++ii) {
            if (xarray[ii] === undefined || yarray[ii] === undefined) continue;
            localCustomLayout.push([
                (xarray[ii]), // * (Pose2D.ZOOM_SPACINGS[0] / Pose2D.ZOOM_SPACINGS[this._zoomLevel]),
                (yarray[ii])// * (Pose2D.ZOOM_SPACINGS[0] / Pose2D.ZOOM_SPACINGS[this._zoomLevel])
            ]);
        }

        // Calculate new base positions
        for (let ii = 0; ii < this._bases.length; ++ii) {
            if (!this._pairs.isPaired(ii)) continue;
            // First work with smaller value, which is either smaller or
            // bigger than the center bp
            this._bases[ii].setXY(
                Math.round(this._bases[ii].x / gridSpace) * gridSpace,
                Math.round(this._bases[ii].y / gridSpace) * gridSpace
            );
            this._bases[ii].setDirty();
        }
        for (let ii = 0; ii < this._bases.length; ++ii) {
            if (!this._pairs.isPaired(ii)) continue;
            for (let jj = 0; jj < localCustomLayout.length; ++jj) {
                localCustomLayout[jj] = [
                    localCustomLayout[jj][0] as number + this._bases[jj].x - this._bases[ii].x,
                    localCustomLayout[jj][1] as number + this._bases[jj].y - this._bases[ii].y
                ];
            }
        }
        // Use setter to force redraw
        this.customLayout = localCustomLayout;
        this._baseRope.enabled = true;
        this._baseRope.redraw(true);
        this._pseudoknotLines.redraw(true);
        this._customLayoutChanged = true;
    }

    public onPoseMouseDown(e: InteractionEvent, closestIndex: number): void {
        const altDown: boolean = Flashbang.app.isAltKeyDown;
        const shiftDown: boolean = Flashbang.app.isShiftKeyDown;
        const ctrlDown: boolean = Flashbang.app.isControlKeyDown || Flashbang.app.isMetaKeyDown;

        if (this.movingAnnotation) {
            return;
        }

        // ctrl + shift: drag base around; ctrl: base mark; shift: shift highlight
        if (closestIndex >= 0) {
            this._mouseDownAltKey = altDown;
            if (ctrlDown && shiftDown) {
                const dragger = new Dragger();
                this.addObject(dragger);

                if (this._currentArrangementTool === Layout.MOVE) {
                    dragger.dragged.connect((p) => {
                        this.onMouseMoved(p as Point, closestIndex);
                    });
                } else if (this._currentArrangementTool === Layout.ROTATE_STEM) {
                    this.rotateStem(closestIndex);
                } else if (this._currentArrangementTool === Layout.FLIP_STEM) {
                    this.flipStem(closestIndex);
                }
                dragger.dragComplete.connect(() => {
                    this.onMouseUp();
                });
                return;
            }
            if (
                (ctrlDown || this.currentColor === RNAPaint.BASE_MARK)
                && closestIndex < this.fullSequenceLength
                && !this.annotationModeActive) {
                this.toggleBaseMark(closestIndex);
                return;
            }
            if (shiftDown && !this.annotationModeActive) {
                if (closestIndex < this.sequenceLength) {
                    this._shiftStart = closestIndex;
                    this._shiftEnd = closestIndex;
                    this.updateShiftHighlight();

                    let reg: Registration | null = null;
                    reg = this.pointerUp.connect(() => {
                        this._shiftStart = -1;
                        this._shiftEnd = -1;
                        if (reg) reg.close();
                    });
                }
                e.stopPropagation();
                return;
            } else if (this.annotationModeActive) {
                if (closestIndex < this.sequenceLength) {
                    let clickedHighlight = false;
                    for (const range of this._annotationRanges) {
                        if (closestIndex >= range.start && closestIndex <= range.end) {
                            clickedHighlight = true;
                            break;
                        } else if (closestIndex <= range.start && closestIndex >= range.end) {
                            clickedHighlight = true;
                            break;
                        }
                    }

                    if (!clickedHighlight) {
                        this._editingAnnotation = true;

                        this._annotationRanges.push({
                            start: closestIndex,
                            end: closestIndex
                        });

                        this.updateAnnotationRangeHighlight();

                        let reg: Registration | null = null;
                        reg = this.pointerUp.connect(() => {
                            this._editingAnnotation = false;
                            if (reg) reg.close();
                        });
                    } else {
                        // Informs others to open dialog
                        this.onCreateAnnotation.value = {
                            ranges: this._annotationRanges
                        };
                    }
                }
                e.stopPropagation();
                return;
            }
            this._lastShiftedCommand = -1;
            this._lastShiftedIndex = -1;
            const cmd: [string, PuzzleEditOp, number[]?] | null = this.parseCommand(this._currentColor, closestIndex);
            if (cmd == null) {
                const dragger = new Dragger();
                this.addObject(dragger);
                dragger.dragged.connect((p) => {
                    this.onMouseMoved(p as Point);
                });
                dragger.dragComplete.connect(() => this.onMouseUp());

                this.onBaseMouseDown(closestIndex, ctrlDown);
            } else {
                this._lastShiftedCommand = this._currentColor;
                this._lastShiftedIndex = closestIndex;

                this.callAddBaseCallback(cmd[0], cmd[1], closestIndex);
            }

            e.stopPropagation();
        } else if (shiftDown && !this.annotationModeActive) {
            this._shiftStart = -1;
            this._shiftEnd = -1;
            this.updateShiftHighlight();
        } else if (this.annotationModeActive) {
            this._annotationRanges = [];
            this._editingAnnotation = false;
            this.updateAnnotationRangeHighlight();
        }
    }

    public toggleBaseMark(baseIndex: number): void {
        if (!this.isTrackedIndex(baseIndex)) {
            this.addBaseMark(baseIndex);
        } else {
            this.removeBaseMark(baseIndex);
        }
    }

    public addBaseMark(baseIndex: number, colors: number | number[] = 0x000000): void {
        if (!this.isTrackedIndex(baseIndex)) {
            if (typeof (colors) === 'number') colors = [colors];
            ROPWait.notifyBlackMark(baseIndex, true);
            this._bases[baseIndex].mark(colors);
        }
    }

    public removeBaseMark(baseIndex: number): void {
        this._bases[baseIndex].unmark();
        ROPWait.notifyBlackMark(baseIndex, false);
    }

    public isTrackedIndex(index: number): boolean {
        return this._bases[index].isMarked();
    }

    public onMouseMoved(point: Point, startIdx?: number): void {
        if (!this._poseField.containsPoint(point.x, point.y)) {
            this.onMouseOut();
            return;
        }

        if (!this._coloring) {
            this.clearMouse();
        }

        this.container.toLocal(point, undefined, Pose2D.P);
        const mouseX = Pose2D.P.x;
        const mouseY = Pose2D.P.y;

        // First, handle the case where you have supplied startIdx, indicating
        // that you are dragging a base to a new location.
        if (startIdx !== undefined) {
            // if (this.customLayout === undefined) {
            const rnaCoords: RNALayout = new RNALayout(
                Pose2D.ZOOM_SPACINGS[this._zoomLevel], Pose2D.ZOOM_SPACINGS[this._zoomLevel]
            );
            rnaCoords.setupTree(this._pairs, this._targetPairs);
            rnaCoords.drawTree(this._customLayout);
            const xarray: number[] = new Array(this._bases.length);
            const yarray: number[] = new Array(this._bases.length);
            rnaCoords.getCoords(xarray, yarray);
            // The simplest thing to do is to use the x/y coords as the new customLayout.
            // This minimizes the calculations you have to do later.
            const localCustomLayout: ([number, number] | [null, null])[] = [];
            for (let ii = 0; ii < this._bases.length; ++ii) {
                if (xarray[ii] === undefined || yarray[ii] === undefined) continue;
                localCustomLayout.push([
                    (xarray[ii]), // * (Pose2D.ZOOM_SPACINGS[0] / Pose2D.ZOOM_SPACINGS[this._zoomLevel]),
                    (yarray[ii])// * (Pose2D.ZOOM_SPACINGS[0] / Pose2D.ZOOM_SPACINGS[this._zoomLevel])
                ]);
            }

            // Ooh, you should drag a helix as a unit.
            if (!this._targetPairs.isPaired(startIdx)) {
                // Update individual base coordinates.
                this._bases[startIdx].setXY(
                    (mouseX - this._offX),
                    (mouseY - this._offY)
                );

                // Update the customLayout in the same way.
                // Actually, after writing this, I no longer know why it works.
                for (let ii = 0; ii < localCustomLayout.length; ++ii) {
                    localCustomLayout[startIdx] = [
                        localCustomLayout[ii][0] as number + (mouseX - this._offX) - this._bases[ii].x,
                        localCustomLayout[ii][1] as number + (mouseY - this._offY) - this._bases[ii].y
                    ];
                }

                this._bases[startIdx].setDirty();
            } else {
                // Find each nt in helix and apply same offset.
                const stem = this._targetPairs.stemWith(startIdx);
                const origX = this._bases[startIdx].x;
                const origY = this._bases[startIdx].y;
                for (const bp of stem) {
                    for (const idx of bp) {
                        this._bases[idx].setXY(
                            mouseX + this._bases[idx].x - origX - this._offX,
                            mouseY + this._bases[idx].y - origY - this._offY
                        );
                        this._bases[idx].setDirty();

                        for (let ii = 0; ii < localCustomLayout.length; ++ii) {
                            localCustomLayout[idx] = [
                                localCustomLayout[ii][0] as number + this._bases[idx].x - this._bases[ii].x,
                                localCustomLayout[ii][1] as number + this._bases[idx].y - this._bases[ii].y
                            ];
                        }
                    }
                }
            }
            // Use setter to force redraw
            this.customLayout = localCustomLayout;
            this._baseRope.enabled = true;
            this._baseRope.redraw(true);
            this._pseudoknotLines.redraw(true);
            this._customLayoutChanged = true;
            return;
        }

        this._paintCursor.display.x = mouseX;
        this._paintCursor.display.y = mouseY;

        let closestDist = -1;
        let closestIndex = -1;
        for (let ii = 0; ii < this.fullSequenceLength; ii++) {
            const mouseDist: number = this._bases[ii].isClicked(
                mouseX - this._offX, mouseY - this._offY, this._zoomLevel, this._coloring
            );
            if (mouseDist >= 0) {
                if (closestIndex < 0 || mouseDist < closestDist) {
                    closestIndex = ii;
                    closestDist = mouseDist;
                }
            }
        }

        if (closestIndex >= 0 && this._currentColor >= 0) {
            this.onBaseMouseMove(closestIndex);
            // document.getElementById(Eterna.PIXI_CONTAINER_ID).style.cursor = 'none';

            if (!this._annotationModeActive) {
                this._paintCursor.display.visible = true;
                this._paintCursor.setShape(this._currentColor);
            }

            const strandName: string | null = this.getStrandName(closestIndex);
            if (strandName != null) {
                this._strandLabel.setText(strandName);
                if (mouseX + 16 + this._strandLabel.width > this._width) {
                    this._strandLabel.display.position = new Point(
                        mouseX - 16 - this._strandLabel.width,
                        mouseY + 16
                    );
                } else {
                    this._strandLabel.display.position = new Point(mouseX + 16, mouseY + 16);
                }
                this._strandLabel.display.visible = true;
            }
        } else {
            this._lastColoredIndex = -1;
        }

        if (!this._coloring) {
            this.updateScoreNodeGui();
        }
    }

    public onMouseUp(): void {
        this.doneColoring();
        this._mouseDownAltKey = false;
        ROPWait.notifyEndPaint();
    }

    public deleteBaseWithIndexPairs(index: number, pairs: SecStruct): [string, PuzzleEditOp, RNABase[]?] {
        if (this.isTrackedIndex(index)) {
            this.toggleBaseMark(index);
        }

        return PoseUtil.deleteNopairWithIndex(index, pairs);
    }

    public clearTracking(): void {
        for (const base of this._bases) {
            base.unmark();
        }
    }

    public get trackedIndices(): { baseIndex: number; colors: number[] }[] {
        const result = [] as { baseIndex: number; colors: number[] }[];
        this._bases.forEach((base, baseIndex) => {
            if (base.isMarked()) {
                result.push({baseIndex, colors: base.markerColors});
            }
        });
        return result;
    }

    public getBase(ind: number): Base {
        return this._bases[ind];
    }

    public get xOffset(): number {
        return this._offX;
    }

    public get yOffset(): number {
        return this._offY;
    }

    public setOffset(offX: number, offY: number): void {
        this._offX = offX;
        this._offY = offY;
        this._redraw = true;
    }

    public get shiftLimit(): number {
        return this._shiftLimit;
    }

    public set shiftLimit(limit: number) {
        this._shiftLimit = limit + this._sequence.length;
    }

    public set barcodes(barcodes: number[]) {
        this._barcodes = barcodes.slice();
    }

    public set puzzleLocks(puzlocks: boolean[] | undefined) {
        if (puzlocks === undefined) {
            this._locks = undefined;
        } else {
            this._locks = puzlocks.slice();
        }

        this._redraw = true;
    }

    public get puzzleLocks(): boolean[] | undefined {
        if (this._locks === undefined) {
            this._locks = Pose2D.createDefaultLocks(this._sequence.length);
        }
        return this._locks.slice();
    }

    public isLocked(seqnum: number): boolean {
        if (this._oligo != null && this._oligoMode === Pose2D.OLIGO_MODE_EXT5P) {
            seqnum -= this._oligo.length;
        }

        if (seqnum < 0 || seqnum >= this._sequence.length) {
            return true;
        } else {
            return this._locks != null && this._locks.length > seqnum && this._locks[seqnum];
        }
    }

    public set forcedStruct(forced: number[] | null) {
        const len: number = this.fullSequenceLength;
        if (forced == null) {
            this._forcedStruct = null;
        } else {
            if (forced.length !== len) {
                throw new Error(`Forced structure length does not match sequence length ${forced.length} ${this._sequence.length} ${this._pairs.length}`);
            }

            this._forcedStruct = forced.slice();
        }

        for (let ii = 0; ii < len; ii++) {
            this._bases[ii].forced = this._forcedStruct != null && this._forcedStruct[ii] !== EPars.FORCE_IGNORE;
        }
    }

    public get forcedStruct(): number[] | null {
        if (this._forcedStruct != null) {
            return this._forcedStruct.slice();
        }

        const temp: number[] = [];
        for (let ii = 0; ii < this.fullSequenceLength; ii++) {
            temp.push(EPars.FORCE_IGNORE);
        }

        return temp;
    }

    public get pseudoknotPairs(): SecStruct {
        return this._pseudoknotPairs;
    }

    public set forcedHighlights(elems: number[]) {
        this._forcedHighlightBox.clear();
        this._forcedHighlightBox.setHighlight(elems);
    }

    public set structConstraints(doCare: boolean[] | undefined) {
        let ii: number;
        const len: number = this.fullSequenceLength;
        const dc: boolean[] | null = (doCare == null ? null : doCare.slice());
        if (dc != null && this._oligosOrder != null) {
            const idxMap: number[] | null = this.getOrderMap(undefined);
            if (idxMap !== null && dc !== null && doCare !== undefined) {
                for (ii = 0; ii < len; ii++) {
                    dc[ii] = doCare[idxMap.indexOf(ii)];
                }
            }
        }
        for (ii = 0; ii < len; ii++) {
            this._bases[ii].dontcare = dc == null ? false : !dc[ii];
        }
    }

    public clearDesignStruct(): void {
        this._designStruct.fill(false);
        this.updateDesignHighlight();
    }

    public toggleDesignStruct(seqnum: number): boolean {
        if (this._designStruct.length !== this.fullSequenceLength) {
            this._designStruct = new Array(this.fullSequenceLength);
        }

        this._designStruct[seqnum] = !this._designStruct[seqnum];
        ROPWait.notifyBlueMark(seqnum, this._designStruct[seqnum]);
        this.updateDesignHighlight();
        const segments: number[] = this.designSegments;
        return (segments.length === 4
            && segments[1] - segments[0] === segments[3] - segments[2]
            && (segments[2] - segments[1] > 3
                || this.fullSequence.hasCut(segments[1], segments[2])));
    }

    public get designSegments(): number[] {
        const elems: number[] = [];
        let curr = 0;
        for (let jj = 0; jj < this.fullSequenceLength; jj++) {
            const stat: number = this._designStruct[jj] ? 1 : 0;
            if ((curr ^ stat) !== 0) {
                elems.push(jj - curr);
                curr = stat;
            }
        }
        if ((elems.length % 2) === 1) {
            elems.push(this.fullSequenceLength - 1);
        }

        return elems;
    }

    public shift3Prime(): void {
        const q: number[] | null = this._shiftHighlightBox.getQueue();
        if (q == null) {
            return;
        }

        const first: number = q[0];
        const last: number = q[1];
        let ii: number;
        // can't shift locked bases
        for (ii = first; ii <= last; ii++) {
            if (this._locks && this._locks[ii]) {
                return;
            }
        }
        // find the next acceptable spot
        let offset = 1;
        const len: number = this.sequenceLength;
        while (last + offset < len) {
            for (ii = first + offset; ii <= last + offset; ii++) {
                if (this._locks && this._locks[ii]) {
                    break;
                }
            }
            if (ii > last + offset) {
                break;
            }
            offset++;
        }
        // if not found, give up
        if (last + offset >= len) {
            return;
        }

        let mutated: RNABase[];
        let segment: RNABase[];
        if (offset === 1) {
            // obtain the segment you are trying to move, plus one 3' base
            segment = this._sequence.baseArray.slice(first, last + 1 + 1);
            // remove the base from the 3' end
            const base = segment.pop();
            Assert.assertIsDefined(base);
            // put the base on the 5' end
            segment.unshift(base);
            mutated = this._sequence.baseArray.slice(0, first)
                .concat(segment)
                .concat(this._sequence.baseArray.slice(last + 1 + 1));
        } else {
            mutated = this._sequence.baseArray.slice();
            for (ii = first; ii <= last; ii++) {
                const xx: number = mutated[ii + offset];
                mutated[ii + offset] = mutated[ii];
                mutated[ii] = xx;
            }
        }

        this._mutatedSequence = this.fullSequence.slice(0);
        this.setMutated(new Sequence(mutated));
        this.doneColoring();
        this._customLayoutChanged = false;
        this._shiftHighlightBox.clear();
        this._shiftHighlightBox.setHighlight([first + offset, last + offset]);
    }

    public shift5Prime(): void {
        const q: number[] | null = this._shiftHighlightBox.getQueue();
        if (q == null) {
            return;
        }

        const first: number = q[0];
        const last: number = q[1];
        let ii: number;
        // can't shift locked bases
        for (ii = first; ii <= last; ii++) {
            if (this._locks && this._locks[ii]) {
                return;
            }
        }
        // find the next acceptable spot
        let ofs = -1;
        while (first + ofs >= 0) {
            for (ii = first + ofs; ii <= last + ofs; ii++) {
                if (this._locks && this._locks[ii]) {
                    break;
                }
            }
            if (ii > last + ofs) {
                break;
            }
            ofs--;
        }
        // if not found, give up
        if (first + ofs < 0) {
            return;
        }

        let mutated: RNABase[];
        let segment: RNABase[];
        if (ofs === -1) {
            segment = this._sequence.baseArray.slice(first - 1, last + 1);
            const base = segment.shift();
            Assert.assertIsDefined(base);
            segment.push(base);
            mutated = this._sequence.baseArray.slice(0, first - 1)
                .concat(segment)
                .concat(this._sequence.baseArray.slice(last + 1));
        } else {
            mutated = this._sequence.baseArray.slice();
            for (ii = first; ii <= last; ii++) {
                const xx: number = mutated[ii + ofs];
                mutated[ii + ofs] = mutated[ii];
                mutated[ii] = xx;
            }
        }

        this._mutatedSequence = this.fullSequence.slice(0);
        this.setMutated(new Sequence(mutated));
        this.doneColoring();
        this._shiftHighlightBox.clear();
        this._shiftHighlightBox.setHighlight([first + ofs, last + ofs]);
    }

    public isDesignStructureHighlighted(index: number): boolean {
        return (this._designStruct[index] === true);
    }

    public getSequenceString(): string {
        return this._sequence.sequenceString();
    }

    public get satisfied(): boolean {
        for (let ii = 0; ii < this._pairs.length; ii++) {
            if (this._pairs.pairingPartner(ii) > ii && !this.isPairSatisfied(ii, this._pairs.pairingPartner(ii))) {
                return false;
            }
        }

        return true;
    }

    public set showNumbering(show: boolean) {
        // FIXME: change  _numberingMode to _showNumbering?
        this._numberingMode = show;
        this._redraw = true;
    }

    public get showNumbering(): boolean {
        return this._numberingMode;
    }

    public set showRope(show: boolean) {
        this._showBaseRope = show;
        this._redraw = true;
    }

    public get showRope(): boolean {
        return this._showBaseRope;
    }

    public set showPseudoknots(show: boolean) {
        this._showPseudoknots = show;
        this._redraw = true;
    }

    public get showPseudoknots(): boolean {
        return this._showPseudoknots;
    }

    public set useSimpleGraphics(simpleGraphics: boolean) {
        this._simpleGraphicsMods = simpleGraphics;
        this._redraw = true;
    }

    public get useSimpleGraphics(): boolean {
        return this._simpleGraphicsMods;
    }

    public set highlightRestricted(highlight: boolean) {
        this._highlightRestricted = highlight;
        this._restrictedHighlightBox.enabled = highlight;
    }

    public get highlightRestricted(): boolean {
        return this._highlightRestricted;
    }

    public set useContinuousExpColors(cont: boolean) {
        this._expContinuous = cont;
        this._redraw = true;

        if (this._expPainter) {
            this.paintFeedback();
        }
    }

    public set useExtendedScale(extended: boolean) {
        this._expExtendedScale = extended;
        this._redraw = true;

        if (this._expPainter) {
            this.paintFeedback();
        }
    }

    public get displayScoreTexts(): boolean {
        return this._displayScoreTexts;
    }

    public set displayScoreTexts(dis: boolean) {
        this._displayScoreTexts = dis;
        this.generateScoreNodes();
    }

    public updateHighlightsAndScores(): void {
        this._prevOffsetX = -1;
        this._prevOffsetY = -1;
        this.generateScoreNodes();
    }

    public set showEnergyHighlight(display: boolean) {
        this._highlightEnergyText = display;
        this.generateScoreNodes();
    }

    public clearHighlight(): void {
        this._selectionHighlightBox.clear();
    }

    // / For restricted queue
    public clearRestrictedHighlight(): void {
        this._restrictedHighlightBox.clear();
    }

    public highlightRestrictedSequence(restricted: number[]): void {
        this._restrictedHighlightBox.setHighlight(restricted);
    }

    public clearUnstableHighlight(): void {
        this._unstableHighlightBox.clear();
    }

    public highlightUnstableSequence(unstable: number[]): void {
        this._unstableHighlightBox.setHighlight(unstable);
    }

    public clearUserDefinedHighlight(): void {
        this._userDefinedHighlightBox.clear();
    }

    public highlightUserDefinedSequence(userDefined: number[]): void {
        this._userDefinedHighlightBox.setHighlight(userDefined);
    }

    public clearShiftHighlight(): void {
        this._shiftHighlightBox.clear();
    }

    public clearAnnotationRanges(): void {
        this._annotationRanges = [];

        this.clearAnnotationHighlight();
    }

    public clearAnnotationHighlight(): void {
        this._annotationHighlightBox.clear();

        if (this._annotationModeActive) {
            for (const base of this._bases) {
                base.container.alpha = Pose2D.ANNOTATION_UNHIGHLIGHTED_OPACITY;
            }
        }
    }

    public praiseStack(stackStart: number, stackEnd: number): void {
        this._praiseQueue.push(stackStart);
        this._praiseQueue.push(stackEnd);
    }

    public praiseSequence(seqStart: number, seqEnd: number): void {
        this._praiseSeq.push(seqStart);
        this._praiseSeq.push(seqEnd);
    }

    private onPraiseStack(stackStart: number, stackEnd: number, playSound: boolean): void {
        let xPos = 0;
        let yPos = 0;

        let playUA = false;
        let playGC = false;
        let playGU = false;

        for (let kk: number = stackStart; kk <= stackEnd; kk++) {
            if (!this._pairs.isPaired(kk)) {
                return;
            }
        }

        for (let ii: number = stackStart; ii <= stackEnd; ii++) {
            const aa: number = ii;
            const bb: number = this._pairs.pairingPartner(ii);

            if ((this._sequence.nt(aa) === RNABase.ADENINE
                && this._sequence.nt(bb) === RNABase.URACIL)
                || (this._sequence.nt(bb) === RNABase.ADENINE
                    && this._sequence.nt(aa) === RNABase.URACIL)) {
                playUA = true;
            } else if ((this._sequence.nt(aa) === RNABase.GUANINE
                && this._sequence.nt(bb) === RNABase.CYTOSINE)
                || (this._sequence.nt(bb) === RNABase.GUANINE
                    && this._sequence.nt(aa) === RNABase.CYTOSINE)) {
                playGC = true;
            } else if ((this._sequence.nt(aa) === RNABase.GUANINE
                && this._sequence.nt(bb) === RNABase.URACIL)
                || (this._sequence.nt(bb) === RNABase.GUANINE
                    && this._sequence.nt(aa) === RNABase.URACIL)) {
                playGU = true;
            }

            this._bases[ii].startSparking();
            this._bases[this._pairs.pairingPartner(ii)].startSparking();
            const p: Point = this.getBaseLoc(ii);
            const p2: Point = this.getBaseLoc(this._pairs.pairingPartner(ii));

            xPos += p.x;
            yPos += p.y;

            xPos += p2.x;
            yPos += p2.y;
        }

        const stackLen: number = (stackEnd - stackStart) + 1;

        xPos /= stackLen * 2;
        yPos /= stackLen * 2;

        const praiseText = stackLen > 1 ? 'Great Pairings!' : 'Great Pairing!';
        const praiseObj = new SceneObject(Fonts.std(praiseText, 20).bold().color(0xffffff).build());
        praiseObj.display.position = new Point(xPos - DisplayUtil.width(praiseObj.display) * 0.5, yPos);
        this.addObject(praiseObj, this.container);

        praiseObj.display.alpha = 0;
        praiseObj.addObject(new SerialTask(
            new ParallelTask(
                new AlphaTask(0.85, 0.33, Easing.easeOut),
                new LocationTask(praiseObj.display.x, praiseObj.display.y - 80, 0.33, Easing.easeOut)
            ),
            new DelayTask(1),
            new ParallelTask(
                new AlphaTask(0, 0.33, Easing.easeOut),
                new LocationTask(praiseObj.display.x, praiseObj.display.y - 120, 0.33, Easing.easeOut)
            ),
            new SelfDestructTask()
        ));

        if (playSound) {
            if (playGC) {
                Flashbang.sound.playSound(Sounds.SoundRG);
            } else if (playUA) {
                Flashbang.sound.playSound(Sounds.SoundYB);
            } else if (playGU) {
                Flashbang.sound.playSound(Sounds.SoundRB);
            }
        }
    }

    public createNewHighlight(nucleotides: number[]): RNAHighlightState {
        const hl: RNAHighlightState = new RNAHighlightState();

        // If any of the nucleotides are part of a stack, highlight its pair as well.
        const addition: number[] = [];
        for (const nuc of nucleotides) {
            if (this._pairs.isPaired(nuc)) {
                addition.push(this._pairs.pairingPartner(nuc));
            }
        }
        nucleotides = nucleotides.concat(addition);

        hl.nuc = nucleotides;
        this._allNewHighlights.push(hl);

        this._redraw = true;
        return hl;
    }

    public removeNewHighlight(highlight: RNAHighlightState): void {
        const idx: number = this._allNewHighlights.indexOf(highlight);
        if (idx >= 0) {
            this._allNewHighlights.splice(idx, 1);
            this._redraw = true;
        }
    }

    private onPraiseSeq(seqStart: number, seqEnd: number): void {
        for (let ii: number = seqStart; ii <= seqEnd; ii++) {
            if (ii >= 0 && ii < this.fullSequenceLength) {
                this._bases[ii].startSparking();
            }
        }
    }

    public startExplosion(): Promise<void> {
        this._isExploding = true;
        this._explosionStartTime = -1;

        if (this._explosionRays.length >= this._sequence.length) {
            for (let ii = 0; ii < this._sequence.length; ii++) {
                const ray = this._explosionRays[ii];
                ray.display.visible = false;
                ray.draw(
                    Vector2.fromPolar(Math.max(this._width, this._height), Math.random() * 2 * Math.PI),
                    this._sequence.nt(ii)
                );
            }
        } else {
            const diff: number = (this._sequence.length - this._explosionRays.length)
                / this._explosionRays.length;
            let diffWalker = 0;
            let rayWalker = 0;

            for (let ii = 0; ii < this._sequence.length; ii++) {
                if (diffWalker < 1) {
                    if (rayWalker >= this._explosionRays.length) {
                        continue;
                    }

                    const ray = this._explosionRays[rayWalker];
                    ray.display.visible = false;
                    ray.draw(
                        Vector2.fromPolar(Math.max(this._width, this._height), Math.random() * 2 * Math.PI),
                        this._sequence.nt(ii)
                    );

                    rayWalker++;
                    diffWalker += diff;
                } else {
                    diffWalker -= 1;
                }
            }
        }

        // If there was an explosion in progress, ensure its promise gets resolved
        if (this._onExplosionComplete != null) {
            this.callExplosionCompleteCallback();
        }

        return new Promise((resolve) => { this._onExplosionComplete = resolve; });
    }

    private callExplosionCompleteCallback(): void {
        if (this._onExplosionComplete != null) {
            const onComplete = this._onExplosionComplete;
            this._onExplosionComplete = null;
            onComplete();
        }
    }

    public clearExplosion(): void {
        if (!this._isExploding) {
            return;
        }

        this._isExploding = false;
        this._explosionStartTime = -1;

        for (const ray of this._explosionRays) {
            ray.fadeOutAndHide();
        }

        this.callExplosionCompleteCallback();
    }

    public set poseEditCallback(cb: () => void) {
        this._poseEditCallback = cb;
    }

    public callPoseEditCallback(): void {
        if (this._poseEditCallback != null) {
            this._poseEditCallback();
        }
    }

    public set trackMovesCallback(cb: (count: number, moves: Move[]) => void) {
        this._trackMovesCallback = cb;
    }

    public callTrackMovesCallback(count: number, moves: Move[]): void {
        if (this._trackMovesCallback != null) {
            this._trackMovesCallback(count, moves);
        }
    }

    public set addBaseCallback(cb: (parenthesis: string | null, op: PuzzleEditOp | null, index: number) => void) {
        this._addBaseCallback = cb;
    }

    public callAddBaseCallback(
        parenthesis: string | null = null, op: PuzzleEditOp | null = null, index: number = -1
    ): void {
        if (this._addBaseCallback != null) {
            this._addBaseCallback(parenthesis, op, index);
        }
    }

    public set startMousedownCallback(cb: PoseMouseDownCallback) {
        this._startMousedownCallback = cb;
    }

    public callStartMousedownCallback(e: InteractionEvent): void {
        e.data.getLocalPosition(this.display, Pose2D.P);
        const mouseX: number = Pose2D.P.x;
        const mouseY: number = Pose2D.P.y;

        let closestDist = -1;
        let closestIndex = -1;

        if (this._startMousedownCallback != null) {
            for (let ii = 0; ii < this.fullSequenceLength; ii++) {
                const mouseDist: number = this._bases[ii].isClicked(
                    mouseX - this._offX, mouseY - this._offY, this._zoomLevel, false
                );
                if (mouseDist >= 0) {
                    if (closestIndex < 0 || mouseDist < closestDist) {
                        closestIndex = ii;
                        closestDist = mouseDist;
                    }
                }
            }
            this._startMousedownCallback(e, closestDist, closestIndex);
        } else {
            this.onPoseMouseDown(e, closestIndex);
        }
    }

    public get satisfiedPairs(): SecStruct {
        return this._pairs.getSatisfiedPairs(this.fullSequence.slice(0));
    }

    public set molecularBindingBonus(bonus: number) {
        this._molecularBindingBonus = bonus;
    }

    public set molecularStructure(pairs: SecStruct | null) {
        if (pairs != null) {
            this._moleculeTargetPairs = pairs.slice(0);
        } else {
            this._moleculeTargetPairs = null;
        }
    }

    public get molecularStructure(): SecStruct | null {
        return this._moleculeTargetPairs;
    }

    public set molecularBindingSite(bindingSite: boolean[] | null) {
        if (bindingSite != null) {
            this._bindingSite = bindingSite.slice();
        } else {
            this._bindingSite = null;
            this.setMolecularBinding(undefined, undefined, this._molecularBindingBonus);
            return;
        }

        const targetPairs: SecStruct | null = this._moleculeTargetPairs
            ? this._moleculeTargetPairs.slice(0)
            : null;
        if (!targetPairs) {
            throw new Error("Can't find molecular target structure");
        }

        const bindingBases: number[] = [];
        const bindingPairs: number[] = [];
        for (let ii = 0; ii < bindingSite.length; ii++) {
            if (bindingSite[ii]) {
                bindingBases.push(ii);
                bindingPairs.push(targetPairs.pairingPartner(ii));
            }
        }
        this.setMolecularBinding(bindingBases, bindingPairs, this._molecularBindingBonus);
    }

    public get molecularBindingSite(): boolean[] | null {
        if (this._bindingSite) {
            return this._bindingSite.slice();
        }

        const temp: boolean[] = [];
        for (let ii = 0; ii < this._sequence.length; ii++) {
            temp.push(false);
        }
        return temp;
    }

    public setMolecularBinding(
        bindingSites: number[] | undefined, bindingPairs: number[] | undefined, bindingBonus: number | undefined
    ): void {
        if (this._molecule != null) {
            this._molecule.destroy({children: true});
            this._molecule = null;
        }

        if (this._molecularBindingBases != null) {
            for (const glow of this._molecularBindingBases) {
                if (glow != null) {
                    glow.destroy({children: true});
                }
            }
            this._molecularBindingBases = null;
        }

        if (bindingSites === undefined || bindingSites.length === 0) {
            return;
        }

        this._molecularBindingBases = new Array(this._sequence.length);
        this._molecularBindingPairs = new Array(this._sequence.length);
        this._molecularBindingBonus = bindingBonus;

        this._molecule = new Molecule();
        this._moleculeLayer.addChild(this._molecule);

        if (bindingPairs === undefined) {
            return;
        }
        for (let ii = 0; ii < bindingSites.length; ii++) {
            const idx = bindingSites[ii];
            const baseGlow = new BaseGlow();
            this._moleculeLayer.addChild(baseGlow);

            this._molecularBindingBases[idx] = baseGlow;
            this._molecularBindingPairs[idx] = bindingPairs[ii];
        }

        this.updateMolecule();
    }

    private updateMolecule(): void {
        if (
            this._molecularBindingBases == null
            || this._molecule == null
        ) {
            return;
        }

        let boundRender = true;
        let boundReal = true;
        const satisfiedPairs: SecStruct = this.satisfiedPairs;

        for (let ii = 0; ii < this._molecularBindingPairs.length; ii++) {
            if (this._molecularBindingBases[ii] == null) {
                continue;
            }
            if (this._molecularBindingPairs[ii] !== this._pairs.pairingPartner(ii)) {
                boundRender = false;
            }

            if (this._molecularBindingPairs[ii] !== satisfiedPairs.pairingPartner(ii)) {
                boundReal = false;
                this._molecularBindingBases[ii].isWrong = true;
            } else {
                this._molecularBindingBases[ii].isWrong = false;
            }
        }
        this._molecule.isWrong = !boundReal;
        this._moleculeIsBound = boundRender;
        this._moleculeIsBoundReal = boundReal;
    }

    public setOligos(oligos?: Oligo[], order?: number[], numPaired: number = 0): void {
        if (oligos === undefined) {
            this._oligos = undefined;
            this._oligosOrder = undefined;
            this._oligosPaired = 0;
            return;
        }

        let same: boolean = this._oligos !== undefined && oligos.length === this._oligos.length;

        if (same) {
            Assert.assertIsDefined(this._oligos);
            for (let k = 0; k < oligos.length && same; k++) {
                if (!Arrays.shallowEqual(this._oligos[k].sequence, oligos[k].sequence)) {
                    same = false;
                    break;
                }
            }
        }

        const prevOrder: number[] | undefined = this._oligosOrder;
        this._oligos = JSON.parse(JSON.stringify(oligos));
        if (order == null) {
            this._oligosOrder = [];
            if (!this._oligos) {
                throw new Error('this._oligos null when we need it not to be!');
            }
            for (let k = 0; k < this._oligos.length; k++) {
                this._oligosOrder[k] = k;
            }
        } else {
            this._oligosOrder = order.slice();
        }
        this._oligosPaired = numPaired;

        const seq: Sequence = this.fullSequence;
        if (seq.length > this._bases.length) {
            const diff: number = (seq.length - this._bases.length);
            for (let k = 0; k < diff; k++) {
                this.createBase();
            }
        }

        const n: number = seq.length;
        for (let k = 0; k < n; k++) {
            this._bases[k].setType(seq.nt(k));
            this._bases[k].baseIndex = k;
        }

        // if possible, maintain visual consistency
        // (strands "fly" from their previous location in the previous oligo order)
        if (same && JSON.stringify(prevOrder) !== JSON.stringify(this._oligosOrder)) {
            const oldX: number[] = [];
            const oldY: number[] = [];
            const idxMap: number[] | null = this.getOrderMap(prevOrder);
            if (idxMap === null) {
                throw new Error('idxMap is null!');
            }
            for (let k = 0; k < seq.length; k++) {
                oldX[k] = this._bases[k].x;
                oldY[k] = this._bases[k].y;
            }
            for (let k = 0; k < seq.length; k++) {
                this._bases[idxMap[k]].setXY(oldX[k], oldY[k]);
            }
        }
    }

    public getOligos(): Oligo[] | null {
        return (this._oligos !== undefined ? JSON.parse(JSON.stringify(this._oligos)) : null);
    }

    public getOrderMap(otherOrder: number[] | undefined): number[] | null {
        if (this._oligos === undefined || this._oligosOrder === undefined) {
            return null;
        }

        const idxMap: number[] = [];
        const ofs: number[] = [];
        let ii: number = this._sequence.length;
        let jj: number;
        for (jj = 0; jj < this._oligos.length; jj++) {
            ofs[this._oligosOrder[jj]] = ii;
            ii += 1 + this._oligos[this._oligosOrder[jj]].sequence.length;
        }
        for (ii = 0; ii < this._sequence.length; ii++) idxMap[ii] = ii;
        for (jj = 0; jj < this._oligos.length; jj++) {
            const zz: number = (otherOrder === undefined ? jj : otherOrder[jj]);
            const kk: number = ofs[zz];
            let xx: number;
            for (xx = 0; xx <= this._oligos[zz].sequence.length; xx++) {
                idxMap[ii + xx] = kk + xx;
            }
            ii += xx;
        }
        return idxMap;
    }

    public saveMarkersContext(): void {
        if (this._oligos === undefined) {
            this._prevOligosOrder = undefined;
        } else if (this._prevOligosOrder === undefined && this._oligosOrder !== undefined) {
            this._prevOligosOrder = this._oligosOrder.slice();
        }
    }

    public transformMarkers(): void {
        if (
            this._prevOligosOrder == null
            || this._oligosOrder == null
            || this._prevOligosOrder.length !== this._oligosOrder.length
        ) {
            this._prevOligosOrder = undefined;
            return;
        }

        const idxMap: number[] | null = this.getOrderMap(this._prevOligosOrder);
        if (idxMap === null) {
            throw new Error('idxMap is null!');
        }
        this._prevOligosOrder = undefined;

        // base marks
        const indices = this.trackedIndices;
        this.clearTracking();
        for (const index of indices) {
            index.baseIndex = idxMap[index.baseIndex];
            this.addBaseMark(index.baseIndex, index.colors);
        }

        // blue highlights ("magic glue")
        const newDesign: boolean[] = [];
        for (let ii = 0; ii < this.fullSequenceLength; ii++) {
            newDesign[idxMap[ii]] = this._designStruct[ii];
        }
        this._designStruct = newDesign;
        this.updateDesignHighlight();
    }

    public setOligo(
        oligo: number[] | undefined,
        mode: number | string | null = Pose2D.OLIGO_MODE_DIMER,
        oName: string | null = null
    ): void {
        if (oligo == null) {
            this._oligo = null;
            return;
        }

        this._oligo = oligo.slice();
        this._oligoName = oName;

        // Puzzle JSON encodes oligoMode as a string, for some reason
        this._oligoMode = typeof (mode) === 'number' ? mode : Number(mode);

        const seq: Sequence = this.fullSequence;
        if (seq.length > this._bases.length) {
            const diff: number = (seq.length - this._bases.length);
            for (let i = 0; i < diff; i++) {
                this.createBase();
            }
        }

        const n: number = seq.length;
        for (let k = 0; k < n; k++) {
            this._bases[k].setType(seq.nt(k));
            this._bases[k].baseIndex = k;
        }
    }

    public set oligoMalus(malus: number) {
        this._oligoMalus = malus;
    }

    public set duplexCost(cost: number) {
        this._duplexCost = cost;
    }

    public set oligoPaired(paired: boolean) {
        const changed: boolean = (this._oligoPaired !== paired);
        this._oligoPaired = paired;
        if (changed) this.updateScoreNodeGui();
    }

    public get fullSequence(): Sequence {
        if (this._oligo == null && this._oligos === undefined) {
            return this._sequence.slice(0);
        }
        let seq: RNABase[] = this._sequence.baseArray.slice();
        if (this._oligos === undefined || this._oligosOrder === undefined) {
            Assert.assertIsDefined(this._oligo);
            if (this._oligoMode === Pose2D.OLIGO_MODE_EXT5P) {
                seq = this._oligo.concat(seq);
            } else {
                if (this._oligoMode === Pose2D.OLIGO_MODE_DIMER) seq.push(RNABase.CUT);
                seq = seq.concat(this._oligo);
            }
            return new Sequence(seq);
        }
        // _oligos != null, we have a multistrand target
        for (let ii = 0; ii < this._oligos.length; ii++) {
            seq.push(RNABase.CUT);
            seq = seq.concat(this._oligos[this._oligosOrder[ii]].sequence);
        }
        return new Sequence(seq);
    }

    public get fullSequenceLength(): number {
        let len: number = this._sequence.length;
        if (this._oligo == null && this._oligos == null) {
            return len;
        }
        if (this._oligos == null && this._oligo !== null) {
            len += this._oligo.length;
            if (this._oligoMode === Pose2D.OLIGO_MODE_DIMER) len++;
            return len;
        }
        Assert.assertIsDefined(this._oligos);
        for (const oligo of this._oligos) {
            len += 1 + oligo.sequence.length;
        }
        return len;
    }

    public getStrandName(seqnum: number): string | null {
        if (this._oligos != null && this._oligosOrder != null && seqnum >= this._sequence.length) {
            let seq: RNABase[] = this._sequence.baseArray.slice();
            for (let ii = 0; ii < this._oligos.length; ii++) {
                seq.push(RNABase.CUT);
                seq = seq.concat(this._oligos[this._oligosOrder[ii]].sequence);
                if (seqnum < seq.length) {
                    let oName: string | undefined = this._oligos[this._oligosOrder[ii]]['name'];
                    if (oName === undefined) oName = `Oligo ${(this._oligosOrder[ii] + 1).toString()}`;
                    return oName;
                }
            }
        }
        if (this._oligo != null && seqnum >= this._sequence.length) {
            return this._oligoName;
        }
        return null;
    }

    public getBoundSequence(): RNABase[] {
        if (this._oligos === undefined || this._oligosOrder === undefined) {
            return this._sequence.baseArray;
        }
        let seq: RNABase[] = this._sequence.baseArray.slice();
        for (let ii = 0; ii < this._oligosPaired; ii++) {
            seq.push(RNABase.CUT);
            seq = seq.concat(this._oligos[this._oligosOrder[ii]].sequence);
        }
        return seq;
    }

    public isPairSatisfied(a: number, b: number): boolean {
        // AMW TODO why swap? do we assume asymmetrical pairs
        if (b < a) {
            const temp: number = a;
            a = b;
            b = temp;
        }

        if (this._pairs.pairingPartner(a) !== b) {
            return false;
        }

        return (EPars.pairType(this.fullSequence.nt(a), this.fullSequence.nt(b)) !== 0);
    }

    public get sequenceLength(): number {
        return this._sequence.length;
    }

    public set sequence(sequence: Sequence) {
        if (Arrays.shallowEqual(this._sequence.baseArray, sequence.baseArray)) {
            return;
        }

        if (this._locks == null) {
            this._locks = Pose2D.createDefaultLocks(this._sequence.length);
        }

        this._sequence = sequence;
        if (this._sequence.length > this._bases.length) {
            const diff: number = (this._sequence.length - this._bases.length);
            for (let ii = 0; ii < diff; ii++) {
                this.createBase();
                this._locks.push(false);
            }
        } else if (this._sequence.length < this._bases.length) {
            for (let ii: number = this._sequence.length; ii < this._bases.length; ii++) {
                if (this.isTrackedIndex(ii)) {
                    this.removeBaseMark(ii);
                }
            }
        }

        const n: number = this.fullSequenceLength;
        for (let ii = 0; ii < n; ii++) {
            if (ii < this._sequence.length) {
                this._bases[ii].setType(this._sequence.nt(ii));
            }
            this._bases[ii].baseIndex = ii;
        }

        this.checkPairs();
        this.updateMolecule();
        this.generateScoreNodes();

        // Update Annotations
        if (this._annotations.length > 0) {
            this.updateAnnotationSpaceAvailability();
            this.eraseAnnotations(true);
            this.drawAnnotations();
        }
    }

    public get sequence(): Sequence {
        return this._sequence.slice(0);
    }

    public getSequenceAt(seq: number): RNABase {
        return this._sequence.nt(seq);
    }

    public set secstruct(pairs: SecStruct) {
        const seq: Sequence = this.fullSequence;
        if (pairs.length !== seq.length) {
            log.debug(pairs.length, seq.length);
            throw new Error("Pair length doesn't match sequence length");
        }

        if (EPars.arePairsSame(pairs, this._pairs)) {
            return;
        }

        this._pairs = pairs.slice(0);

        // AMW: We don't have to worry about this case where... pairs are
        // asymmetric somehow?
        // for (let ii = 0; ii < this._pairs.length; ii++) {
        //     if (this._pairs.pairingPartner(ii) > ii) {
        //         this._pairs.pairingPartner(this._pairs.pairingPartner(ii) = ii;
        //     }
        // }

        // Recompute sequence layout
        this.computeLayout(false);
        this.checkPairs();
        this.updateMolecule();
        this.generateScoreNodes();

        // Update Annotations
        if (this._annotations.length > 0) {
            this.updateAnnotationSpaceAvailability();
            this.eraseAnnotations(true);
            this.drawAnnotations();
        }
    }

    public get secstruct(): SecStruct {
        return this._pairs.slice(0);
    }

    public set targetPairs(setting: SecStruct) {
        this._targetPairs = setting.slice(0);
        // for (let ii = 0; ii < this._targetPairs.length; ii++) {
        //     // AMW TODO: symmetrizing the secstruct directly; seems not ideal.
        //     if (this._targetPairs.pairingPartner(ii) > ii) {
        //         this._targetPairs.pairingPartner(this._targetPairs.pairingPartner(ii)) = ii;
        //     }
        // }
    }

    public set customLayout(setting: Array<[number, number] | [null, null]> | undefined) {
        // Compare first
        if (setting === undefined) {
            if (this.customLayout !== undefined) {
                this._customLayout = setting;
                this.computeLayout(false);
            }
        } else if (this.customLayout === undefined) {
            this._customLayout = setting;
            this.computeLayout(false);
        } else {
            let diff = false;
            for (let ii = 0; ii < setting.length; ++ii) {
                if (setting[ii] !== this._customLayout?.[ii] ?? null) {
                    diff = true;
                    break;
                }
            }
            if (!diff) return;
            this._customLayout = setting;
            this.computeLayout(false);
        }
    }

    public get customLayout(): Array<[number, number] | [null, null]> | undefined {
        return this._customLayout;
    }

    public set customNumbering(setting: (number | null)[] | undefined) {
        this._customNumbering = setting;
    }

    public get customNumbering(): (number | null)[] | undefined {
        return this._customNumbering;
    }

    public set pseudoknotted(pk: boolean) {
        this._pseudoknotted = pk;
    }

    public get pseudoknotted(): boolean {
        return this._pseudoknotted;
    }

    public checkOverlap(): boolean {
        const radius: number = Pose2D.ZOOM_SPACINGS[0];
        const rnaDrawer: RNALayout = new RNALayout(radius, radius);
        rnaDrawer.setupTree(this._pairs, this._targetPairs);
        rnaDrawer.drawTree(this._customLayout);

        const xarray: number[] = new Array(this._bases.length);
        const yarray: number[] = new Array(this._bases.length);
        rnaDrawer.getCoords(xarray, yarray);
        for (let ii = 0; ii < this._bases.length; ii++) {
            const ax: number = xarray[ii];
            const ay: number = yarray[ii];
            for (let jj: number = ii + 2; jj < this._bases.length; jj++) {
                let bx: number = xarray[jj];
                let by: number = yarray[jj];
                bx = ax - bx;
                by = ay - by;
                if (bx * bx + by * by < (radius * radius) / 10) {
                    return true;
                }
            }
        }
        return false;
    }

    // highlight the base before the cursor
    public trackCursor(index: number | null): void {
        this._cursorIndex = index;
        if (this._cursorIndex !== null && this._cursorIndex > 0) {
            const center: Point = this.getBaseLoc(this._cursorIndex - 1);
            if (this._cursorBox == null) {
                this._cursorBox = new Graphics();
                this.container.addChild(this._cursorBox);
            }
            this._cursorBox.x = center.x;
            this._cursorBox.y = center.y;
            this._cursorBox.visible = true;
            this._cursorBox.clear();
            this._cursorBox.lineStyle(Base.MARKER_THICKNESS * Base.MARKER_RADIUS[this.zoomLevel],
                Pose2D.COLOR_CURSOR);
            this._cursorBox.drawCircle(0, 0, Base.MARKER_RADIUS[this.zoomLevel]);
        } else if (this._cursorBox != null) {
            this._cursorBox.destroy({children: true});
            this._cursorBox = null;
        }
    }

    public get trackedCursorIdx(): number | null {
        return this._cursorIndex;
    }

    private makeHighlightState() {
        let hlState: RNAHighlightState | undefined;
        if (this._allNewHighlights.length > 0) {
            hlState = new RNAHighlightState();
            hlState.nuc = [];
            hlState.isOn = true;
            for (const existingHighlight of this._allNewHighlights) {
                if (existingHighlight.nuc !== null) {
                    hlState.nuc = hlState.nuc.concat(existingHighlight.nuc);
                }
            }
        }
        return hlState;
    }

    private setAllDrawParams(fullSeq: Sequence, currentTime: number, hlState?: RNAHighlightState): void {
        for (let ii = 0; ii < fullSeq.length; ii++) {
            // skip the oligo separator
            if (fullSeq.nt(ii) === RNABase.CUT) {
                continue;
            }

            const useBarcode = (this._barcodes != null && this._barcodes.indexOf(ii) >= 0);

            this._bases[ii].forceUnpaired = (
                this._forcedStruct != null && this._forcedStruct[ii] === EPars.FORCE_UNPAIRED
            );

            const drawFlags: number = BaseDrawFlags.builder()
                .locked(this.isLocked(ii))
                .letterMode(this._lettermode)
                .lowPerform(this._simpleGraphicsMods)
                .useBarcode(useBarcode)
                .result();

            let numberBitmap: Texture | null = null;
            if (this._numberingMode) {
                let displayNumber: number | null = ii + 1;
                if (this._customNumbering != null) displayNumber = this._customNumbering[ii];
                if ((displayNumber != null)
                    && (ii === 0 || displayNumber % 5 === 0 || ii === fullSeq.length - 1)) {
                    numberBitmap = BitmapManager.getNumberBitmap(displayNumber);
                }
            }

            this._bases[ii].setDrawParams(
                this._zoomLevel, this._offX, this._offY, currentTime, drawFlags, numberBitmap, hlState
            );
        }
    }

    /* override */
    public update(_dt: number): void {
        if (!this.display.worldVisible) {
            // update is expensive, so don't bother doing it if we're not visible
            return;
        }
        Assert.assertIsDefined(this.mode);
        const currentTime: number = this.mode.time;
        for (const anchor of this._anchoredObjects) {
            if (anchor.isLive) {
                const p: Point = this.getBaseLoc(anchor.base);
                anchor.object.display.position = new Point(p.x + anchor.offset.x, p.y + anchor.offset.y);
            }
        }

        const fullSeq: Sequence = this.fullSequence;
        let center: Point;

        // Hide bases that aren't part of our current sequence
        if (!this._showNucleotideRange) {
            for (let ii = 0; ii < this._bases.length; ++ii) {
                this._bases[ii].display.visible = this.isNucleotidePartOfSequence(ii);
            }
        }

        const basesMoved = this._baseToX && this._baseToY && this._baseFromX && this._baseFromY;
        if (basesMoved) {
            // Update base locations

            if (this._foldStartTime < 0) {
                this._foldStartTime = currentTime;
            }

            let prog = (currentTime - this._foldStartTime) / (this._foldDuration);

            if (prog >= 1) {
                prog = 1;
                this._offsetTranslating = false;

                // Update Annotations
                if (this._annotations.length > 0) {
                    this.updateAnnotationSpaceAvailability();
                    this.eraseAnnotations(true);
                    this.drawAnnotations();
                }
            }

            if (this._offsetTranslating) {
                this._redraw = true;
                this._offX = prog * this._endOffsetX + (1 - prog) * this._startOffsetX;
                this._offY = prog * this._endOffsetY + (1 - prog) * this._startOffsetY;
            }

            this.setAnimationProgress(prog);
        } else if (currentTime - this.lastSampledTime > 2 && !this._isExploding) {
            this.lastSampledTime = currentTime;

            for (let ii = 0; ii < fullSeq.length; ii++) {
                if (!this._pairs.isPaired(ii) && !this._simpleGraphicsMods && Math.random() > 0.7) {
                    this._bases[ii].animate();
                }
            }
        }

        // / Update score node
        if (!this._annotationModeActive) {
            this.updateScoreNodeVisualization(this._offX !== this._prevOffsetX || this._offY !== this._prevOffsetY);
        }

        // / Bitblt rendering
        const needRedraw = this._bases.some(
            (base) => base.needRedraw(this._simpleGraphicsMods)
        );

        if (needRedraw || this._redraw) {
            // Create highlight state to pass to bases, then set up draw params.
            this.setAllDrawParams(fullSeq, currentTime, this.makeHighlightState());
        }

        // AMW TODO: this means that PuzzleEditMode can get a baserope showing
        // if the custom layout tools are used.
        this._baseRope.enabled = this._showBaseRope || (this._customLayout != null);
        this._pseudoknotLines.enabled = this._pseudoknotPairs
            && this._pseudoknotPairs.nonempty();

        if (this._redraw || basesMoved) {
            this._baseRope.redraw(true /* force baseXY */);
            if (this.pseudoknotPairs && this.pseudoknotPairs.length !== 0) {
                this._pseudoknotLines.redraw(true /* force baseXY */);
            }

            if (this._cursorIndex != null && this._cursorIndex > 0 && this._cursorBox !== null) {
                center = this.getBaseLoc(this._cursorIndex - 1);
                this._cursorBox.x = center.x;
                this._cursorBox.y = center.y;
                this._cursorBox.visible = true;
                this._cursorBox.clear();
                this._cursorBox.lineStyle(Base.MARKER_THICKNESS * Base.MARKER_RADIUS[this.zoomLevel],
                    Pose2D.COLOR_CURSOR);
                this._cursorBox.drawCircle(0, 0, Base.MARKER_RADIUS[this.zoomLevel]);
            }
        }

        this._redraw = false;

        this._moleculeLayer.visible = false;
        if (this._molecularBindingBases != null && this._molecule != null) {
            let molX = 0;
            let molY = 0;
            let nbases = 0;
            for (let ii = 0; ii < fullSeq.length; ii++) {
                const baseglow: BaseGlow = this._molecularBindingBases[ii];
                if (baseglow != null) {
                    const pos: Point = this._bases[ii].getLastDrawnPos();
                    baseglow.updateView(this._zoomLevel, pos.x, pos.y, currentTime);
                    molX += pos.x;
                    molY += pos.y;
                    nbases += 1;
                }
            }

            if (nbases > 0) {
                molX /= nbases;
                molY /= nbases;
            }

            if (!this._moleculeIsBound) {
                molX = 30;
                molY = 200;
            }

            this._molecule.updateView(this._zoomLevel, molX, molY, currentTime);
            this._moleculeLayer.visible = true;
        }

        if (fullSeq.findCut() >= 0) {
            if (this._oligoBases == null) {
                this._oligoBases = new Array(fullSeq.length);
            }

            const boundLen: number = this.getBoundSequence().length;
            for (let ii = fullSeq.findCut() + 1; ii < fullSeq.length; ii++) {
                if (this._oligoBases[ii] === undefined) {
                    this._oligoBases[ii] = new BaseGlow();
                }
                const baseglow = this._oligoBases[ii];
                if ((this._oligoPaired || (this._oligosPaired > 0 && ii < boundLen)) && this._pairs.isPaired(ii)) {
                    baseglow.isWrong = this._restrictedHighlightBox.isInQueue(ii);
                    const pos = this._bases[ii].getLastDrawnPos();
                    baseglow.updateView(this._zoomLevel, pos.x, pos.y, currentTime);
                    this._moleculeLayer.visible = true;
                }
            }
        }

        let goX = 0;
        let goY = 0;

        for (let ii = 0; ii < fullSeq.length - 1; ii++) {
            let outX: number = goX;
            let outY: number = goY;
            let dirSign = 1;
            if (ii < this._baseRotationDirectionSign.length) dirSign = this._baseRotationDirectionSign[ii];

            if (this._sequence.length < fullSeq.length && ii === this._sequence.length - 1) {
                this._bases[ii].setGoDir(goX, goY);
                this._bases[ii].setOutDir(dirSign * -goY, dirSign * goX);
                this._bases[ii].setLast(true);
                continue;
            }

            goX = this._bases[ii + 1].x - this._bases[ii].x;
            goY = this._bases[ii + 1].y - this._bases[ii].y;

            const goLength: number = Math.sqrt(goX * goX + goY * goY);
            if (goLength > Pose2D.ZOOM_SPACINGS[this._zoomLevel]) {
                goX *= (Pose2D.ZOOM_SPACINGS[this._zoomLevel] / goLength);
                goY *= (Pose2D.ZOOM_SPACINGS[this._zoomLevel] / goLength);
            }

            outX += goX;
            outY += goY;

            if (ii > 0) {
                outX /= 2.0;
                outY /= 2.0;
            }

            this._bases[ii].setGoDir(goX, goY);
            this._bases[ii].setOutDir(dirSign * -outY, dirSign * outX);
            this._bases[ii].setLast(false);
        }

        if (fullSeq.length >= 1) {
            this._bases[fullSeq.length - 1].setGoDir(goX, goY);
            this._bases[fullSeq.length - 1].setOutDir(-goY, goX);
            this._bases[fullSeq.length - 1].setLast(true);
        }

        // / Praise stacks when RNA is not moving
        if (!this._offsetTranslating && this._baseToX == null) {
            if (this._praiseQueue.length > 0) {
                for (let ii = 0; ii < this._praiseQueue.length; ii += 2) {
                    this.onPraiseStack(
                        this._praiseQueue[ii],
                        this._praiseQueue[ii + 1],
                        (ii + 1) === (this._praiseQueue.length - 1)
                    );
                }
                this._praiseQueue = [];
            } else if (this._praiseSeq.length > 0) {
                for (let ii = 0; ii < this._praiseSeq.length; ii += 2) {
                    this.onPraiseSeq(this._praiseSeq[ii], this._praiseSeq[ii + 1]);
                }
                this._praiseSeq = [];
            }
        }

        if (this._isExploding && !this._offsetTranslating && this._baseToX == null) {
            if (this._explosionStartTime < 0) {
                this._explosionStartTime = currentTime;
                this._origOffsetX = this._offX;
                this._origOffsetY = this._offY;
            }

            this._offX = this._origOffsetX + (Math.random() * 2 - 1) * 5;
            this._offY = this._origOffsetY + (Math.random() * 2 - 1) * 5;
            this._redraw = true;

            const prog = (currentTime - this._explosionStartTime) * 5;

            if (this._explosionRays.length >= fullSeq.length) {
                for (let ii = 0; ii < Math.min(prog, fullSeq.length); ii++) {
                    const ray = this._explosionRays[ii];

                    if (!ray.display.visible) {
                        ray.display.alpha = 0;
                        ray.display.visible = true;
                        ray.fadeIn();
                    }

                    ray.display.position = this.getBaseLoc(ii);
                }
            } else {
                const diff: number = (fullSeq.length - this._explosionRays.length) / this._explosionRays.length;
                let diffWalker = 0;
                let rayWalker = 0;

                for (let ii = 0; ii < fullSeq.length; ii++) {
                    if (diffWalker < 1) {
                        if (rayWalker >= this._explosionRays.length || rayWalker >= prog) {
                            continue;
                        }

                        const ray = this._explosionRays[rayWalker];
                        if (!ray.display.visible) {
                            ray.display.alpha = 0;
                            ray.display.visible = true;
                            ray.fadeIn();
                        }

                        ray.display.position = this.getBaseLoc(ii);
                        rayWalker++;
                        diffWalker += diff;
                    } else {
                        diffWalker -= 1;
                    }
                }
            }

            if (prog >= Math.min(fullSeq.length, this._explosionRays.length) + 10) {
                this.clearExplosion();
            }
        }

        this._prevOffsetX = this._offX;
        this._prevOffsetY = this._offY;
    }

    public setAnimationProgress(progress: number) {
        if (this._baseToX && this._baseToY && this._baseFromX && this._baseFromY) {
            for (let ii = 0; ii < this.fullSequence.length; ii++) {
                const vx: number = this._baseToX[ii] - this._baseFromX[ii];
                const vy: number = this._baseToY[ii] - this._baseFromY[ii];

                const currentX: number = this._baseFromX[ii] + ((vx + (vx * progress)) / 2) * progress;
                const currentY: number = this._baseFromY[ii] + ((vy + (vy * progress)) / 2) * progress;

                this._bases[ii].setXY(currentX, currentY);
            }
        }

        if (progress >= 1) {
            this._baseToX = null;
            this._baseToY = null;
            this._baseFromX = null;
            this._baseFromY = null;

            this.updateScoreNodeGui();

            if (this.checkOverlap()) {
                // If overlaps have been introduced, make sure the explosion factor input is shown
                this._explosionFactorPanel.display.visible = true;
            } else if (this._explosionFactorPanel.display.visible === true) {
                // If all overlaps have been removed, remove the explosion
                this._explosionFactor = 1;
                this._explosionFactorPanel.display.visible = false;
                this.computeLayout(true);
                this._redraw = true;
            }
        }
    }

    public numPairs(satisfied: boolean): number {
        // AMW TODO: this is very similar to SecStruct::numPairs, but with satisfied.
        let n = 0;
        for (let ii = 0; ii < this._pairs.length; ii++) {
            if (this._pairs.pairingPartner(ii) > ii
                    && (!satisfied || this.isPairSatisfied(ii, this._pairs.pairingPartner(ii)))) {
                n++;
            }
        }
        return n;
    }

    public set lettermode(lettermode: boolean) {
        this._lettermode = lettermode;
        this._redraw = true;
    }

    public get lettermode(): boolean {
        return this._lettermode;
    }

    public get showTotalEnergy(): boolean {
        return this._showTotalEnergy;
    }

    public set showTotalEnergy(show: boolean) {
        this._showTotalEnergy = show;
        this._primaryScoreEnergyDisplay.visible = (show && this._scoreFolder != null);
        this._secondaryScoreEnergyDisplay.visible = (
            show && this._scoreFolder != null && this._secondaryScoreEnergyDisplay.hasText
        );
        this._deltaScoreEnergyDisplay.visible = show && this._scoreFolder != null;
    }

    public get showExplosionFactor(): boolean {
        return this._explosionFactorPanel.display.visible;
    }

    public set showExplosionFactor(show: boolean) {
        this._explosionFactorPanel.display.visible = show;
    }

    public set scoreFolder(folder: Folder | null) {
        if (this._scoreFolder !== folder) {
            this._scoreFolder = folder;
            this.showTotalEnergy = this._showTotalEnergy;
            this.generateScoreNodes();
        }
    }

    public baseShiftWithCommand(command: number, index: number): void {
        const cmd: [string, PuzzleEditOp, number[]?] | null = this.parseCommand(command, index);
        if (cmd != null) {
            const parenthesis: string = cmd[0];
            const op: PuzzleEditOp = cmd[1];
            this.baseShift(parenthesis, op, index);
        }
    }

    public baseShift(parenthesis: string, op: PuzzleEditOp, index: number): void {
        let sequence: RNABase[] = this.sequence.baseArray;
        let locks: boolean[] | undefined = this.puzzleLocks;
        let bindingSite: boolean[] | null = this.molecularBindingSite;
        const sequenceBackup: RNABase[] = this.sequence.baseArray;
        const locksBackup: boolean[] | undefined = this.puzzleLocks;
        const bindingSiteBackup: boolean[] | null = this.molecularBindingSite;
        let pindex: number;

        if (sequence.length > parenthesis.length) {
            sequence = sequence.slice(0, parenthesis.length);
            locks = locks ? locks.slice(0, parenthesis.length) : undefined;
            bindingSite = bindingSite ? bindingSite.slice(0, parenthesis.length) : null;
        }

        for (let ii: number = sequence.length; ii < parenthesis.length; ii++) {
            sequence.push(RNABase.ADENINE);
            if (locks) locks.push(false);
            if (bindingSite) bindingSite.push(false);
        }
        // BASE SHIFTING MODIFIED HERE. Delete comments to apply the changes
        if (op === PuzzleEditOp.ADD_BASE) {
            // Add a base
            const afterIndex: number[] = sequence.slice(index);
            const afterLockIndex: boolean[] | null = locks ? locks.slice(index) : null;
            const afterBindingSiteIndex: boolean[] | null = bindingSite ? bindingSite.slice(index) : null;

            sequence[index] = RNABase.ADENINE;
            if (locks) locks[index] = false;
            if (bindingSite) bindingSite[index] = false;

            for (let ii = 0; ii < afterIndex.length - 1; ii++) {
                sequence[ii + index + 1] = afterIndex[ii];
                if (locks && afterLockIndex) locks[ii + index + 1] = afterLockIndex[ii];
                if (bindingSite && afterBindingSiteIndex) bindingSite[ii + index + 1] = afterBindingSiteIndex[ii];
            }
        } else if (op === PuzzleEditOp.ADD_PAIR) {
            // Add a pair
            pindex = this.secstruct.pairingPartner(index);
            const afterIndex = sequence.slice(index);
            const afterLockIndex = locks ? locks.slice(index) : null;
            const afterBindingSiteIndex = bindingSite ? bindingSite.slice(index) : null;

            sequence[index] = RNABase.ADENINE;
            sequence[pindex + 2] = RNABase.ADENINE;
            if (locks) locks[index] = false;
            if (locks) locks[pindex + 2] = false;
            if (bindingSite) bindingSite[index] = false;
            if (bindingSite) bindingSite[pindex + 2] = false;

            for (let ii = 0; ii < afterIndex.length - 2; ii++) {
                if (ii + index > pindex) {
                    sequence[ii + index + 2] = afterIndex[ii];
                    if (locks && afterLockIndex) locks[ii + index + 2] = afterLockIndex[ii];
                    if (bindingSite && afterBindingSiteIndex) bindingSite[ii + index + 2] = afterBindingSiteIndex[ii];
                } else {
                    sequence[ii + index + 1] = afterIndex[ii];
                    if (locks && afterLockIndex) locks[ii + index + 1] = afterLockIndex[ii];
                    if (bindingSite && afterBindingSiteIndex) bindingSite[ii + index + 1] = afterBindingSiteIndex[ii];
                }
            }
        } else if (op === PuzzleEditOp.ADD_CYCLE) {
            // Add a cycle of length 3
            const afterIndex = sequence.slice(index);
            const afterLockIndex = locks ? locks.slice(index) : null;
            const afterBindingSiteIndex = bindingSite ? bindingSite.slice(index) : null;

            sequence[index] = RNABase.ADENINE;
            sequence[index + 1] = RNABase.ADENINE;
            sequence[index + 2] = RNABase.ADENINE;
            sequence[index + 3] = RNABase.ADENINE;
            sequence[index + 4] = RNABase.ADENINE;

            if (locks) {
                locks[index] = false;
                locks[index + 1] = false;
                locks[index + 2] = false;
                locks[index + 3] = false;
                locks[index + 4] = false;
            }

            if (bindingSite) {
                bindingSite[index] = false;
                bindingSite[index + 1] = false;
                bindingSite[index + 2] = false;
                bindingSite[index + 3] = false;
                bindingSite[index + 4] = false;
            }

            for (let ii = 0; ii < afterIndex.length - 5; ii++) {
                sequence[ii + index + 5] = afterIndex[ii];
                if (locks && afterLockIndex) locks[ii + index + 5] = afterLockIndex[ii];
                if (bindingSite && afterBindingSiteIndex) bindingSite[ii + index + 5] = afterBindingSiteIndex[ii];
            }
        } else if (op === PuzzleEditOp.DELETE_PAIR) {
            // Delete a pair
            pindex = this.secstruct.pairingPartner(index);
            const afterIndex = sequenceBackup.slice(index + 1);
            const afterLockIndex = locksBackup ? locksBackup.slice(index + 1) : null;
            const afterBindingSiteIndex = bindingSiteBackup ? bindingSiteBackup.slice(index + 1) : null;

            for (let ii = 0; ii < afterIndex.length - 1; ii++) {
                if (ii + index >= pindex - 1) {
                    sequence[ii + index] = afterIndex[ii + 1];
                    if (locks && afterLockIndex) locks[ii + index] = afterLockIndex[ii + 1];
                    if (bindingSite && afterBindingSiteIndex) bindingSite[ii + index] = afterBindingSiteIndex[ii + 1];
                } else {
                    sequence[ii + index] = afterIndex[ii];
                    if (locks && afterLockIndex) locks[ii + index] = afterLockIndex[ii];
                    if (bindingSite && afterBindingSiteIndex) bindingSite[ii + index] = afterBindingSiteIndex[ii];
                }
            }
        } else if (op === PuzzleEditOp.DELETE_BASE) {
            // Delete a base
            const afterIndex = sequenceBackup.slice(index + 1);
            const afterLockIndex = locksBackup ? locksBackup.slice(index + 1) : null;
            const afterBindingSiteIndex = bindingSiteBackup ? bindingSiteBackup.slice(index + 1) : null;

            for (let ii = 0; ii < afterIndex.length; ii++) {
                sequence[ii + index] = afterIndex[ii];
                if (locks && afterLockIndex) locks[ii + index] = afterLockIndex[ii];
                if (bindingSite && afterBindingSiteIndex) bindingSite[ii + index] = afterBindingSiteIndex[ii];
            }
        }

        this.sequence = new Sequence(sequence);
        this.puzzleLocks = locks;
        this.molecularStructure = SecStruct.fromParens(parenthesis);
        this.molecularBindingSite = bindingSite;
    }

    public registerPaintTool(paintColor: number, tool: Booster): void {
        this._dynPaintColors.push(paintColor);
        this._dynPaintTools.push(tool);
    }

    public get lastShiftedIndex(): number {
        return this._lastShiftedIndex;
    }

    public get lastShiftedCommand(): number {
        return this._lastShiftedCommand;
    }

    public setBaseColor(seqpos: number, inColor: RNABase): void {
        this._mutatedSequence = this._sequence.slice(0);
        this._mutatedSequence.setNt(seqpos, inColor);
        this._bases[seqpos].setType(inColor, true);

        this._lastColoredIndex = seqpos;
        this._bases[seqpos].animate();
        this.doneColoring();
    }

    public forceEditable(b: boolean, editList: number[] | null = null): void {
        this._editable = b;
        this._editableIndices = editList;
    }

    /**
     * Center a nucleotide into view
     * @param index: 1-based index of the nucleotide, or its custom display number
     *
     */
    public focusNucleotide(index: number) {
        const baseIndex = (() => {
            if (this._customNumbering) {
                return this._customNumbering.findIndex((e) => e === index);
            } else {
                return index - 1;
            }
        })();

        if (baseIndex < 0 || baseIndex >= this._bases.length) {
            // eslint-disable-next-line
            console.warn(`Can't focus nucleotide with index '${index}'`);
            return;
        }

        this.setOffset(
            this._width / 2 - this._bases[baseIndex].x,
            this._height / 2 - this._bases[baseIndex].y
        );
    }

    public showNucleotideRange(range: [number, number] | null) {
        const [start, end] = range ?? [1, this._bases.length];
        if (start < 1 || end > this._bases.length || start >= end) {
            // eslint-disable-next-line
            console.warn(`Invalid nucleotide range [${start}, ${end}]`);
            return;
        }

        this._showNucleotideRange = Boolean(range);
        if (!range) {
            return;
        }

        for (let i = 0; i < this._bases.length; ++i) {
            this._bases[i].container.visible = i >= (start - 1)
                && i < end
                && this.isNucleotidePartOfSequence(i);
        }
    }

    private computeLayout(fast: boolean = false): void {
        const fullSeq: Sequence = this.fullSequence;

        if (fullSeq.length > this._bases.length) {
            log.debug(fullSeq.length, this._bases.length);
            throw new Error("Sequence length and pose length don't match");
        }

        const n: number = fullSeq.length;
        let xMid = 0;
        let yMid = 0;

        let xarray: number[] = new Array(n);
        let yarray: number[] = new Array(n);

        let exceptionIndices: number[] | undefined;
        if (fullSeq.findCut() >= 0) {
            exceptionIndices = [];
            exceptionIndices.push(0);
            let oligoIndex = -1;
            // array of positions of connectors "&"
            while (fullSeq.findCut(oligoIndex + 1) >= 0) {
                oligoIndex = fullSeq.findCut(oligoIndex + 1);
                exceptionIndices.push(oligoIndex);
            }
        }
        const rnaDrawer: RNALayout = new RNALayout(
            Pose2D.ZOOM_SPACINGS[this._zoomLevel],
            Pose2D.ZOOM_SPACINGS[this._zoomLevel] * this._explosionFactor,
            exceptionIndices
        );

        rnaDrawer.setupTree(this._pairs, this._targetPairs);
        rnaDrawer.drawTree(this._customLayout);
        rnaDrawer.getCoords(xarray, yarray);
        this._pseudoknotPairs = rnaDrawer.pseudoknotPairs;

        this._baseRotationDirectionSign = new Array(n);
        rnaDrawer.getRotationDirectionSign(this._baseRotationDirectionSign);

        if (this._desiredAngle === 90) {
            const tmp = xarray;
            xarray = yarray;
            yarray = tmp;
        }

        const xmin: number = Math.min(...xarray);
        const xmax: number = Math.max(...xarray);
        const ymin: number = Math.min(...yarray);
        const ymax: number = Math.max(...yarray);

        xMid = (xmax + xmin) / 2.0;
        yMid = (ymax + ymin) / 2.0;

        this._baseFromX = this._bases.map((b) => b.x);
        this._baseFromY = this._bases.map((b) => b.y);
        this._baseToX = xarray.map((x) => x - xMid);
        this._baseToY = yarray.map((y) => y - yMid);

        this._foldStartTime = -1;
        if (fast) {
            this._foldDuration = 0.45;
        } else {
            this._foldDuration = 0.7;
        }
    }

    private onMouseOut(): void {
        this.clearMouse();
        this.updateScoreNodeGui();
    }

    private deleteBaseWithIndex(index: number): [string, PuzzleEditOp, RNABase[]?] {
        if (this.isTrackedIndex(index)) {
            this.toggleBaseMark(index);
        }

        if (!this._pairs.isPaired(index) || this.isLocked(this._pairs.pairingPartner(index))) {
            return PoseUtil.deleteNopairWithIndex(index, this._pairs);
        } else {
            return PoseUtil.deletePairWithIndex(index, this._pairs);
        }
    }

    private onBaseMouseDown(seqnum: number, togglelock: boolean): void {
        this._lastColoredIndex = seqnum;

        if (togglelock || !this.isEditable(seqnum)) return;

        this._coloring = true;
        this._mutatedSequence = this.fullSequence.slice(0);

        if (this._currentColor === RNAPaint.LOCK) {
            if (!this._locks) {
                this._locks = [];
                for (let ii = 0; ii < this._sequence.length; ii++) {
                    this._locks.push(false);
                }
            }
            this._locks[seqnum] = !this._locks[seqnum];
            this._bases[seqnum].setDirty();
            this._lockUpdated = true;
        } else if (this._currentColor === RNAPaint.BINDING_SITE) {
            if (this._bindingSite != null && this._bindingSite[seqnum]) {
                this._bindingSite = [];
                for (let ii = 0; ii < this._sequence.length; ii++) {
                    this._bindingSite.push(false);
                }
                this.molecularBindingSite = this._bindingSite;
                this._bindingSiteUpdated = true;
            } else {
                const bindingBases: number[] | null = this._pairs.isInternal(seqnum);
                if (bindingBases != null && bindingBases.length > 4) {
                    this._bindingSite = [];
                    for (let ii = 0; ii < this._sequence.length; ii++) {
                        this._bindingSite.push(false);
                    }

                    for (let ii = 0; ii < bindingBases.length; ii++) {
                        this._bindingSite[bindingBases[ii]] = true;
                    }
                    this.molecularBindingSite = this._bindingSite;
                    this._bindingSiteUpdated = true;
                } else {
                    (this.mode as GameMode).showNotification(
                        'Binding site can be only formed at loops between 2 stacks\n(Internal loops and Bulges)'
                    );
                }
            }
        } else if (this._mouseDownAltKey || this._currentColor === RNAPaint.MAGIC_GLUE) {
            if (this.toggleDesignStruct(seqnum)) {
                this._designStructUpdated = true;
            }
        } else if (!this.isLocked(seqnum)) {
            if (this._currentColor >= 1 && this._currentColor <= 4) {
                this._mutatedSequence.setNt(seqnum, this._currentColor);
                ROPWait.notifyPaint(seqnum, this._bases[seqnum].type, this._currentColor);
                this._bases[seqnum].setType(this._currentColor, true);
            } else if (this._currentColor === RNAPaint.PAIR && this._pairs.isPaired(seqnum)) {
                const pi = this._pairs.pairingPartner(seqnum);
                if (this.isLocked(pi)) {
                    return;
                }

                const clickBase: RNABase = this._mutatedSequence.nt(seqnum);

                this._mutatedSequence.setNt(seqnum, this._mutatedSequence.nt(pi));
                this._mutatedSequence.setNt(pi, clickBase);

                this._bases[seqnum].setType(this._mutatedSequence.nt(seqnum), true);
                this._bases[pi].setType(this._mutatedSequence.nt(pi), true);
            } else if (this._currentColor === RNAPaint.AU_PAIR && this._pairs.isPaired(seqnum)) {
                const pi = this._pairs.pairingPartner(seqnum);
                if (this.isLocked(pi)) {
                    return;
                }

                this._mutatedSequence.setNt(seqnum, RNABase.ADENINE);
                this._mutatedSequence.setNt(pi, RNABase.URACIL);

                this._bases[seqnum].setType(this._mutatedSequence.nt(seqnum), true);
                this._bases[pi].setType(this._mutatedSequence.nt(pi), true);
            } else if (this._currentColor === RNAPaint.GC_PAIR && this._pairs.isPaired(seqnum)) {
                const pi = this._pairs.pairingPartner(seqnum);
                if (this.isLocked(pi)) {
                    return;
                }

                this._mutatedSequence.setNt(seqnum, RNABase.GUANINE);
                this._mutatedSequence.setNt(pi, RNABase.CYTOSINE);

                this._bases[seqnum].setType(this._mutatedSequence.nt(seqnum), true);
                this._bases[pi].setType(this._mutatedSequence.nt(pi), true);
            } else if (this._currentColor === RNAPaint.GU_PAIR && this._pairs.isPaired(seqnum)) {
                const pi = this._pairs.pairingPartner(seqnum);
                if (this.isLocked(pi)) {
                    return;
                }

                this._mutatedSequence.setNt(seqnum, RNABase.URACIL);
                this._mutatedSequence.setNt(pi, RNABase.GUANINE);

                this._bases[seqnum].setType(this._mutatedSequence.nt(seqnum), true);
                this._bases[pi].setType(this._mutatedSequence.nt(pi), true);
            } else if (this._dynPaintColors.indexOf(this._currentColor) >= 0) {
                const index: number = this._dynPaintColors.indexOf(this._currentColor);
                this._dynPaintTools[index].onPaint(this, seqnum);
            }
        }
    }

    private onBaseMouseMove(seqnum: number): void {
        if (!this._coloring && this._shiftStart >= 0 && seqnum < this.sequenceLength) {
            this._shiftEnd = seqnum;
            this.updateShiftHighlight();
        } else if (
            !this._coloring
            && this.annotationModeActive
            && this._annotationRanges.length > 0
            && this._editingAnnotation
            && seqnum < this.sequenceLength
        ) {
            this._annotationRanges[this._annotationRanges.length - 1].end = seqnum;
            this.updateAnnotationRangeHighlight();
        }

        if (!this._coloring || (seqnum === this._lastColoredIndex)) {
            return;
        }

        if (this._currentColor === RNAPaint.LOCK) {
            if (!this._locks) {
                this._locks = [];
                for (let ii = 0; ii < this._sequence.length; ii++) {
                    this._locks.push(false);
                }
            }
            this._locks[seqnum] = !this._locks[seqnum];
            this._bases[seqnum].setDirty();
            this._lockUpdated = true;
        } else if (this._mouseDownAltKey || this._currentColor === RNAPaint.MAGIC_GLUE) {
            if (this.toggleDesignStruct(seqnum)) {
                this._designStructUpdated = true;
            }
        } else if (!this.isLocked(seqnum)) {
            if (this._mutatedSequence === null) {
                throw new Error('The clicked base is not locked, but the mutated sequence is null: critical error!');
            }
            if (this._currentColor >= 1 && this._currentColor <= 4) {
                this._mutatedSequence.setNt(seqnum, this._currentColor);
                ROPWait.notifyPaint(seqnum, this._bases[seqnum].type, this._currentColor);
                this._bases[seqnum].setType(this._currentColor, true);
            } else if (this._currentColor === RNAPaint.PAIR) {
                if (this._pairs.isPaired(seqnum)) {
                    const pi = this._pairs.pairingPartner(seqnum);

                    if (this.isLocked(pi)) {
                        return;
                    }

                    const clickBase: number = this._mutatedSequence.nt(seqnum);

                    this._mutatedSequence.setNt(seqnum, this._mutatedSequence.nt(pi));
                    this._mutatedSequence.setNt(pi, clickBase);

                    this._bases[seqnum].setType(this._mutatedSequence.nt(seqnum), true);
                    this._bases[pi].setType(this._mutatedSequence.nt(pi), true);
                }
            } else if (this._currentColor === RNAPaint.AU_PAIR) {
                if (this._pairs.isPaired(seqnum)) {
                    const pi = this._pairs.pairingPartner(seqnum);

                    if (this.isLocked(pi)) {
                        return;
                    }

                    this._mutatedSequence.setNt(seqnum, RNABase.ADENINE);
                    this._mutatedSequence.setNt(pi, RNABase.URACIL);

                    this._bases[seqnum].setType(this._mutatedSequence.nt(seqnum), true);
                    this._bases[pi].setType(this._mutatedSequence.nt(pi), true);
                }
            } else if (this._currentColor === RNAPaint.GC_PAIR) {
                if (this._pairs.isPaired(seqnum)) {
                    const pi = this._pairs.pairingPartner(seqnum);

                    if (this.isLocked(pi)) {
                        return;
                    }

                    this._mutatedSequence.setNt(seqnum, RNABase.GUANINE);
                    this._mutatedSequence.setNt(pi, RNABase.CYTOSINE);

                    this._bases[seqnum].setType(this._mutatedSequence.nt(seqnum), true);
                    this._bases[pi].setType(this._mutatedSequence.nt(pi), true);
                }
            } else if (this._currentColor === RNAPaint.GU_PAIR) {
                if (this._pairs.isPaired(seqnum)) {
                    const pi = this._pairs.pairingPartner(seqnum);

                    if (this.isLocked(pi)) {
                        return;
                    }

                    this._mutatedSequence.setNt(seqnum, RNABase.URACIL);
                    this._mutatedSequence.setNt(pi, RNABase.GUANINE);

                    this._bases[seqnum].setType(this._mutatedSequence.nt(seqnum), true);
                    this._bases[pi].setType(this._mutatedSequence.nt(pi), true);
                }
            } else if (this._dynPaintColors.indexOf(this._currentColor) >= 0) {
                const index: number = this._dynPaintColors.indexOf(this._currentColor);
                this._dynPaintTools[index].onPainting(this, seqnum);
            }
        }
        this._lastColoredIndex = seqnum;
        this._bases[seqnum].animate();
    }

    private updateDesignHighlight(): void {
        const elems: number[] = this.designSegments;
        this._selectionHighlightBox.clear();
        this._selectionHighlightBox.setHighlight(elems);
    }

    private updateShiftHighlight(): void {
        this._shiftHighlightBox.clear();
        if (this._shiftStart >= 0) {
            this._shiftHighlightBox.setHighlight(
                this._shiftEnd < this._shiftStart
                    ? [this._shiftEnd, this._shiftStart] : [this._shiftStart, this._shiftEnd]
            );
        }
    }

    private setAnnotationRangeHighlight(ranges: AnnotationRange[]): void {
        this._annotationHighlightBox.clear();
        for (const range of ranges) {
            this._annotationHighlightBox.setHighlight([range.start, range.end]);
        }
    }

    private updateAnnotationRangeHighlight(): void {
        this._annotationHighlightBox.clear();

        if (this._annotationModeActive) {
            for (let i = 0; i < this._bases.length; i++) {
                this._bases[i].container.alpha = Pose2D.ANNOTATION_UNHIGHLIGHTED_OPACITY;
            }
        }

        if (this._annotationRanges.length > 0) {
            for (const range of this._annotationRanges) {
                this._annotationHighlightBox.setHighlight(
                    range.end < range.start
                        ? [range.end, range.start] : [range.start, range.end]
                );

                if (range.end >= range.start) {
                    for (let i = range.start; i <= range.end; i++) {
                        this._bases[i].container.alpha = 1;
                    }
                } else {
                    for (let i = range.end; i <= range.start; i++) {
                        this._bases[i].container.alpha = 1;
                    }
                }
            }
        }
    }

    public updateAnnotations(annotations: AnnotationData[]): void {
        this.clearAnnotationHighlight();

        const updatedAnnotations = [];
        for (const annotation of annotations) {
            const existingAnnotationIndex = this._annotations.findIndex(
                (obj: AnnotationDisplayObject) => annotation.id === obj.data.id
            );
            if (existingAnnotationIndex !== -1) {
                // Keep existing positions in case we have custom positioning
                updatedAnnotations.push({
                    data: annotation,
                    type: AnnotationItemType.ANNOTATION,
                    positions: this._annotations[existingAnnotationIndex].positions,
                    displays: []
                });
            } else {
                updatedAnnotations.push({
                    data: annotation,
                    type: AnnotationItemType.ANNOTATION,
                    positions: [],
                    displays: []
                });
            }

            if (annotation.selected && annotation.visible && annotation.ranges) {
                this.setAnnotationRangeHighlight(annotation.ranges);
            }
        }

        this._annotations = updatedAnnotations;

        // Update Annotations
        // We don't check for annotations.length > 0 because
        // we want to account for scenario where to go
        // from non-zero to zero annotation
        if (this._annotationSpaceAvailability.length === 0) {
            this.updateAnnotationSpaceAvailability();
        }
        this.eraseAnnotations(true);
        this.drawAnnotations();
    }

    public updateLayers(layers: AnnotationData[]): void {
        this.clearAnnotationHighlight();

        this._layers = [];
        for (const layer of layers) {
            this._layers.push({
                data: layer,
                type: AnnotationItemType.LAYER,
                positions: [],
                displays: []
            });

            let ranges: AnnotationRange[] = [];
            if (layer.children) {
                for (const annotation of layer.children) {
                    if (annotation.ranges) {
                        ranges = ranges.concat(annotation.ranges);
                    }
                }
            }

            if (layer.selected) {
                this.setAnnotationRangeHighlight(ranges);
            }
        }

        // Update Annotations
        // We don't check for annotations.length > 0 because
        // we want to account for scenario where to go
        // from non-zero to zero annotation
        this.eraseAnnotations();
        this.drawAnnotations();
    }

    public updateGraph(collection: AnnotationDataCollection): void {
        const structurePath = (node: AnnotationData): AnnotationGraphNode => {
            let positions: AnnotationPosition[] = [];
            if (node.type === AnnotationItemType.ANNOTATION) {
                // Annotation
                const annotationIndex = this._annotations.findIndex(
                    (obj: AnnotationDisplayObject) => node.id === obj.data.id
                );

                if (annotationIndex !== -1) {
                    positions = this._annotations[annotationIndex].positions;
                }
            } else {
                // Layer
                const annotationIndex = this._layers.findIndex(
                    (obj: AnnotationDisplayObject) => node.id === obj.data.id
                );

                if (annotationIndex !== -1) {
                    positions = this._annotations[annotationIndex].positions;
                }
            }

            const children: AnnotationGraphNode[] = [];
            if (node.children) {
                for (const child of node.children) {
                    children.push(structurePath(child));
                }
            }

            return {
                data: node,
                type: node.type,
                positions,
                children
            };
        };

        const annotationGraph: AnnotationGraph = {
            puzzle: [],
            solution: []
        };
        const puzzleGraph: AnnotationGraphNode[] = [];
        for (const child of collection.puzzle) {
            const path = structurePath(child);
            puzzleGraph.push(path);
        }
        annotationGraph.puzzle = puzzleGraph;
        const solutionGraph: AnnotationGraphNode[] = [];
        for (const child of collection.solution) {
            const path = structurePath(child);
            solutionGraph.push(path);
        }
        annotationGraph.solution = solutionGraph;

        this._annotationGraph = annotationGraph;
    }

    public get annotationGraph() {
        return this._annotationGraph;
    }

    public set puzzleAnnotationsEditable(editable: boolean) {
        this._puzzleAnnotationsEditable = editable;
    }

    public get puzzleAnnotationsEditable(): boolean {
        return this._puzzleAnnotationsEditable;
    }

    public set annotationModeActive(active: boolean) {
        this._annotationModeActive = active;
        this._redraw = true;
        if (!active) {
            this.clearAnnotationHighlight();
            // Change base cursors
            for (const base of this._bases) {
                base.container.alpha = 1;
            }

            // Make annotation canvas opaque
            if (this._annotations.length > 0) {
                this.annotationCanvas.alpha = 1;
            }
        } else {
            // Change base cursors
            for (const base of this._bases) {
                base.container.alpha = Pose2D.ANNOTATION_UNHIGHLIGHTED_OPACITY;
            }

            // Make annotation canvas translucent
            if (this._annotations.length > 0) {
                this.annotationCanvas.alpha = Pose2D.ANNOTATION_UNHIGHLIGHTED_OPACITY;
            }
        }
    }

    public get annotationModeActive(): boolean {
        return this._annotationModeActive;
    }

    public updateAnnotationSpaceAvailability(): void {
        // Set annotation space availability to true
        this._annotationSpaceAvailability = Array(this._baseLayer.height).fill(0).map(
            () => Array(this._baseLayer.width).fill(true)
        );

        const baseLayerBounds = DisplayUtil.getBoundsRelative(this._baseLayer, this.container);

        // Populate with bases
        for (let i = 0; i < this._bases.length; i++) {
            const base = this._bases[i];
            const baseBounds = DisplayUtil.getBoundsRelative(base.display, this._baseLayer);
            const baseRowStart = Math.max(0, Math.floor(baseBounds.y - baseLayerBounds.y));
            const baseRowEnd = Math.min(
                Math.ceil(this._baseLayer.height),
                Math.ceil(baseBounds.y - baseLayerBounds.y + baseBounds.height)
            );
            const baseColStart = Math.max(0, Math.floor(baseBounds.x - baseLayerBounds.x));
            const baseColEnd = Math.min(
                Math.ceil(this._baseLayer.width),
                Math.ceil(baseBounds.x - baseLayerBounds.x + baseBounds.width)
            );

            for (let row = baseRowStart; row < baseRowEnd; row++) {
                const replaceCount = baseColEnd - baseColStart;
                this._annotationSpaceAvailability[row].splice(
                    baseColStart,
                    replaceCount,
                    ...Array(replaceCount).fill(false)
                );
            }
        }
    }

    private static drawEnergyHighlight(hilite: Graphics, energy: Sprite): Graphics {
        // Draw highlight around the energy reading.
        // Give it a bit of padding so the highlight isn't so tight.
        const PADDING = 2;
        return hilite.clear()
            .lineStyle(1, 0xFFFFFF, 0.7)
            .drawRoundedRect(
                energy.x - PADDING, energy.y - PADDING,
                energy.width + PADDING, energy.height + PADDING, 10
            );
    }

    private updateEnergyHighlight(energy: Sprite, idx: number, vis: boolean): void {
        if (idx >= this._energyHighlights.length) {
            if (!this._highlightEnergyText) {
                return;
            }

            const obj = new SceneObject(Pose2D.drawEnergyHighlight(new Graphics(), energy));
            obj.display.alpha = 0;
            obj.addObject(new RepeatingTask((): SerialTask => new SerialTask(
                new AlphaTask(1, 0.5),
                new AlphaTask(0, 0.5)
            )));

            this.addObject(obj, this.container);
            this._energyHighlights.push(obj);
        } else {
            const obj = this._energyHighlights[idx];
            obj.display.visible = vis;
            Pose2D.drawEnergyHighlight(obj.display as Graphics, energy);
        }
    }

    private clearEnergyHighlights(): void {
        for (const obj of this._energyHighlights) {
            obj.destroySelf();
        }

        this._energyHighlights = [];
    }

    public eraseAnnotations(reset: boolean = false, ignoreCustom: boolean = false): void {
        if (this.annotationCanvas.children.length > 0) {
            // Clear prior annotation displays
            this._annotations.forEach((annotation) => {
                annotation.displays.forEach((display) => display.destroySelf);
                annotation.displays = [];
            });

            this._resetAnnotationPositions = reset;
            this._ignoreCustomAnnotationPositions = ignoreCustom;

            // Remove from any remaining artifacts from canvas
            this.annotationCanvas.removeChildren();
            this.annotationCanvas.clear();
        }
    }

    /**
     * Renders annotation cards near ranges of interest
     */
    public drawAnnotations(): void {
        const getAnnotationCard = (item: AnnotationDisplayObject): AnnotationCard => {
            let textColor;
            switch (item.data.category) {
                case AnnotationCategory.STRUCTURE:
                    textColor = AnnotationItem.STRUCTURE_RIBBON_COLOR;
                    break;
                case AnnotationCategory.PUZZLE:
                    textColor = AnnotationItem.PUZZLE_RIBBON_COLOR;
                    break;
                default:
                    textColor = AnnotationItem.SOLUTION_RIBBON_COLOR;
                    break;
            }

            const card = new AnnotationCard(
                item.type,
                item.data,
                this._puzzleAnnotationsEditable,
                textColor
            );

            card.pointerOver.connect(() => {
                // highlight associated range
                if (item.type === AnnotationItemType.ANNOTATION) {
                    const annotation = item.data as AnnotationData;
                    if (annotation.ranges) {
                        this.setAnnotationRangeHighlight(annotation.ranges);
                    }
                } else if (item.type === AnnotationItemType.LAYER) {
                    const layer = item.data as AnnotationData;
                    let ranges: AnnotationRange[] = [];
                    if (layer.children) {
                        for (const annotation of layer.children) {
                            if (annotation.ranges) {
                                ranges = ranges.concat(annotation.ranges);
                            }
                        }
                    }
                    this.setAnnotationRangeHighlight(ranges);
                }
            });
            card.pointerOut.connect(() => {
                // remove associated range
                if (!item.data.selected) {
                    this.clearAnnotationHighlight();
                }
            });

            card.pointerDown.connect(() => {
                if (item.type === AnnotationItemType.ANNOTATION) {
                    this.onSelectAnnotation.value = item.data as AnnotationData;
                } else if (item.type === AnnotationItemType.LAYER) {
                    this.onSelectLayer.value = item.data as AnnotationData;
                }
            });
            card.isMoving.connect((moving: boolean) => {
                this.movingAnnotation = moving;
            });

            if (item.type === AnnotationItemType.ANNOTATION) {
                // We don't need to apply access control logic here
                // This is handled in AnnotationCard
                card.onEditButtonPressed.connect(() => {
                    this.onEditAnnotation.value = item.data as AnnotationData;
                    this.onEditAnnotation.value = null;
                });
            }

            return card;
        };

        /**
         * Attempts to place a single annotation item
         *
         * @param item display object data with positioning and annotation metadata
         * @param itemIndex index of item within parent array
         */
        const placeItem = (item: AnnotationDisplayObject, itemIndex: number): void => {
            // Skip if annotation is marked as hidden
            if (!item.data.visible) return;

            // If annotation positions have been computed already
            // use cached value
            if (
                item.positions.length > 0
                && !this._resetAnnotationPositions
                && !this._ignoreCustomAnnotationPositions
            ) {
                for (let i = 0; i < item.positions.length; i++) {
                    const position = item.positions[i];
                    const card = getAnnotationCard(item);
                    if (item.type === AnnotationItemType.ANNOTATION) {
                        card.onMovedAnnotation.connect((point: Point) => {
                            const anchorIndex = this._annotations[itemIndex].positions[i].anchorIndex;
                            const anchorPoint = new Point(
                                this._bases[anchorIndex].x + this._offX,
                                this._bases[anchorIndex].y + this._offY
                            );
                            // Compute relative position
                            this._annotations[itemIndex].positions[i].relPosition = new Point(
                                point.x - anchorPoint.x,
                                point.y - anchorPoint.y
                            );
                            this._annotations[itemIndex].positions[i].custom = true;
                        });
                    }
                    item.displays.push(card);
                    this.addObject(card, this.annotationCanvas);

                    const anchorIndex = position.anchorIndex;
                    const anchorPoint = new Point(
                        this._bases[anchorIndex].x + this._offX,
                        this._bases[anchorIndex].y + this._offY
                    );
                    card.display.position = new Point(
                        position.relPosition.x + anchorPoint.x,
                        position.relPosition.y + anchorPoint.y
                    );
                }

                return;
            }

            let ranges: AnnotationRange[] = [];
            if (item.type === AnnotationItemType.LAYER) {
                // We only want one layer label for each item, so we pick the first range we find
                //
                // An improvement that could be made is to find the
                // "center of mass" of all the ranges in a layer
                // and place the layer label at an appropriate base closest to
                // the center of mass point
                const layerData = item.data as AnnotationData;

                if (layerData.children) {
                    for (const annotation of layerData.children) {
                        if (annotation.ranges) {
                            ranges.push(annotation.ranges[0]);
                            break;
                        }
                    }
                }
            } else {
                const annotationData = item.data as AnnotationData;
                if (annotationData.ranges) {
                    ranges = annotationData.ranges;
                }
            }

            // Future Improvement:
            // A single annotation or layer can be associated with multiple ranges
            // We handle this by generating a label for each range, regardless
            // of the proximity of the ranges
            //
            // An improvement that can be made is to only "duplicate" labels
            // if range positions exceed some defined threshold to avoid unnecessary
            // label duplicates.
            for (let i = 0; i < ranges.length; i++) {
                const range = ranges[i];
                const prevPosition = item.positions.length > i ? item.positions[i] : null;

                const card = getAnnotationCard(item);
                if (item.type === AnnotationItemType.ANNOTATION) {
                    card.onMovedAnnotation.connect((point: Point) => {
                        const anchorIndex = this._annotations[itemIndex].positions[i].anchorIndex;
                        const anchorPoint = new Point(
                            this._bases[anchorIndex].x + this._offX,
                            this._bases[anchorIndex].y + this._offY
                        );
                        // Compute relative position
                        this._annotations[itemIndex].positions[i].relPosition = new Point(
                            point.x - anchorPoint.x,
                            point.y - anchorPoint.y
                        );
                        this._annotations[itemIndex].positions[i].custom = true;
                    });
                }
                // We need to prematruely add this to the display object graph
                // so that we can read it's dimensions/position
                this.addObject(card, this.annotationCanvas);

                let absolutePosition: Point | null = null;
                let relPosition: Point | null = null;
                let anchorIndex: number | null = null;
                let customPosition = false;
                let zoomLevel: number = this._zoomLevel;
                let preventCustomPositionOverwrite = false;
                if (prevPosition?.custom && !this._ignoreCustomAnnotationPositions) {
                    relPosition = prevPosition.relPosition;
                    anchorIndex = prevPosition.anchorIndex;
                    const zoomScaling = 1
                    + 2 * ((prevPosition.zoomLevel - this._zoomLevel) / Pose2D.ZOOM_SPACINGS.length);
                    const anchorPoint = new Point(
                        this._bases[anchorIndex].x + this._offX,
                        this._bases[anchorIndex].y + this._offY
                    );
                    absolutePosition = new Point(
                        relPosition.x * zoomScaling + anchorPoint.x,
                        relPosition.y * zoomScaling + anchorPoint.y
                    );
                    zoomLevel = prevPosition.zoomLevel;
                    customPosition = true;
                } else {
                    // Make anchor midpoint of range
                    // Account for reverse ranges
                    if (range.start < range.end) {
                        anchorIndex = range.start + Math.floor((range.end - range.start) / 2);
                    } else {
                        anchorIndex = range.end + Math.floor((range.start - range.end) / 2);
                    }

                    // Make sure anchor sits within sequence length
                    if (anchorIndex >= this._bases.length - 1) continue;

                    const anchorPoint = new Point(
                        this._bases[anchorIndex].x + this._offX,
                        this._bases[anchorIndex].y + this._offY
                    );

                    // Run a search to find best place to locate
                    // annotation within available space
                    relPosition = this.computeAnnotationPositionPoint(
                        anchorIndex,
                        anchorIndex,
                        anchorPoint,
                        this._bases[anchorIndex].display,
                        card,
                        0
                    );

                    if (relPosition) {
                        absolutePosition = new Point(
                            relPosition.x + anchorPoint.x,
                            relPosition.y + anchorPoint.y
                        );
                    }

                    if (prevPosition?.custom && this._ignoreCustomAnnotationPositions) {
                        // If we're ignoring custom annotation positions
                        // we still want to differentiate between custom and non-custom
                        // but we don't want to overwrite the position as we would
                        // when not ignoring
                        preventCustomPositionOverwrite = true;
                    }
                }

                // Handle position
                if (relPosition && absolutePosition) {
                    // Set position
                    card.display.position = absolutePosition;
                    // Save display
                    item.displays.push(card);

                    // Cache position
                    if (item.displays.length > item.positions.length) {
                        item.positions.push({
                            anchorIndex,
                            relPosition,
                            zoomLevel,
                            custom: customPosition
                        });
                    } else if (!preventCustomPositionOverwrite) {
                        // We don't want to over
                        item.positions[i] = {
                            anchorIndex,
                            relPosition,
                            zoomLevel,
                            custom: customPosition
                        };
                    }
                } else {
                    // We should ideally always receive a position
                    //
                    // In cases we dont, remove annotation card from view
                    card.destroySelf();
                }
            }
        };

        if (this._zoomLevel > Pose2D.ANNOTATION_LAYER_THRESHOLD) {
            // visualize layers
            for (let i = 0; i < this._layers.length; i++) {
                placeItem(this._layers[i], i);
            }
        } else {
            // visualize annotations
            for (let i = 0; i < this._annotations.length; i++) {
                placeItem(this._annotations[i], i);
            }
        }

        // Make canvas translucent if we have annotation mode active
        if (this._annotationModeActive) {
            this.annotationCanvas.alpha = Pose2D.ANNOTATION_UNHIGHLIGHTED_OPACITY;
        } else {
            this.annotationCanvas.alpha = 1;
        }

        this._resetAnnotationPositions = false;
    }

    /**
     * Searches for spot in available space around annnotation range
     * to place annotation card
     *
     * Runs a recursive/iterative search on each place defined about the co-ordinate
     * system with the anchor point as the origin until it finds a place
     * @param originalAnchorIndex the index of the base associated with the initial call
     * @param currentAnchorIndex the index of the base currently being used as the anchor
     * @param anchorPoint mid-point of anchor that defines the origin on which calculations are made
     * relative from. We use the mid-point and not the top-left corner, as is convention in pixi,
     * because bases in Eterna.js have their position saved as a central point
     * @param anchorDisplay display object of anchor
     * @param annotationCard annotation to be placed
     * @param numSearchAttempts the number of new search attempts undergone
     * @param anchorOffsetX the x-offset applied to the anchor point when attempting a recursive call.
     * The effect is to treat an annotation card as the new anchor while still maintaining reference
     * to the anchor of interest.
     * @param anchorOffsetY the y-offset applied to the anchor point when attempting a recursive call.
     * The effect is to treat an annotation card as the new anchor while still maintaining reference
     * to the anchor of interest.
     * @param includeCenters whether we attempt to place the annotation vertically centered on the
     * left or right of origin
     * @return relative position (to anchor) of annotation relative computed from the annotation's
     * top-left corner
     */
    private computeAnnotationPositionPoint(
        originalAnchorIndex: number,
        currentAnchorIndex: number,
        anchorPoint: Point,
        anchorDisplay: Container,
        annotationCard: AnnotationCard,
        numSearchAttempts: number,
        anchorOffsetX: number = 0,
        anchorOffsetY: number = 0,
        includeCenters: boolean = true
    ): Point | null {
        // DEBUG SETTINGS:
        // In order to visualize placement conflicts, you can toggle these variables
        // RED RECTANGLES = Placement conflicts with bases
        // ORANGE RECTANGLES = Placement conflicts with existing annotations
        // GREEN RECTANGLES = Placement availability
        const DEBUG_ANNOTATION_PLACEMENT = false;
        const DEBUG_FROM_SEARCH_ATTEMPT = 0;
        const BASE_CONFLICT_COLOR = 0xFF0000;
        const ANNOTATION_CONFLICT_COLOR = 0xFF6500;
        const POSSIBLE_PLACEMENT_AVAILABILITY_COLOR = 0x00FF00;

        // x and y offsets that create space between annotations
        const xOffset: number = Pose2D.DEFAULT_ANNOTATION_SHIFT;
        const yOffset = 0;

        /**
         * Determines which place defined about a co-ordinate
         * system with an anchor point as the origin are occupied by bases
         *
         *                   top center
         *     top-left   --------------- top-right
         *               |       |       |
         *               |       |       |
         *               |       |       |
         *   left-center  ---- Anchor ---  right-center
         *               |       |       |
         *               |       |       |
         *               |       |       |
         *   bottom-left  ---------------  bottom-right
         *                 bottom center
         *
         * @return array of occupied places
         */
        const findBaseConflicts = (): AnnotationBaseConflicts => {
            /**
             * Inspects a particular place for placement availability
             *
             * @param startRow pixel row from which to begin vertical inspection
             * @param stopRow pixel row from which to end vertical inspection
             * @param startCol pixel column from which to begin horizontal inspection
             * @param stopCol pixel column from which to end horizontal inspection
             */
            const inspectRegion = (
                startRow: number,
                stopRow: number,
                startCol: number,
                stopCol: number
            ): AnnotationBaseConflict | null => {
                let conflict: AnnotationBaseConflict | null = null;
                for (let row = Math.floor(
                    Math.min(Math.max(0, startRow), this._annotationSpaceAvailability.length)
                );
                    row < Math.ceil(
                        Math.min(Math.max(0, stopRow), this._annotationSpaceAvailability.length)
                    );
                    row++) {
                    for (let col = Math.floor(
                        Math.min(Math.max(0, startCol), this._annotationSpaceAvailability[0].length)
                    );
                        col < Math.ceil(
                            Math.min(Math.max(0, stopCol), this._annotationSpaceAvailability[0].length)
                        );
                        col++) {
                        // build out conflict bounds
                        if (!this._annotationSpaceAvailability[row][col] && !conflict) {
                            conflict = {
                                bounds: new Rectangle(
                                    col,
                                    row,
                                    1,
                                    1
                                ),
                                resolvable: false
                            };
                        } else if (!this._annotationSpaceAvailability[row][col] && conflict) {
                            conflict.bounds = new Rectangle(
                                conflict.bounds.x,
                                conflict.bounds.y,
                                Math.max(conflict.bounds.width, col - conflict.bounds.x),
                                Math.max(conflict.bounds.height, row - conflict.bounds.y)
                            );
                        }
                    }
                }

                return conflict;
            };

            /**
             * Searches for any corrections (minor position shifts) that could be applied
             * to resolve a placement confict
             *
             * @param conflict bounds where conflict occurs
             * @param startRow pixel row from which to begin vertical inspection
             * @param stopRow pixel row from which to end vertical inspection
             * @param startCol pixel column from which to begin horizontal inspection
             * @param stopCol pixel column from which to end horizontal inspection
             * @param testUp whether to search in the upwards direction
             * @param testDown whether to search in the downwards direction
             * @param testLeft whether to search in the leftwards direction
             * @param testRight whether to search in the rightwards direction
             */
            const testForConflictCorrection = (
                conflict: AnnotationBaseConflict,
                startRow: number,
                stopRow: number,
                startCol: number,
                stopCol: number,
                testUp: boolean,
                testDown: boolean,
                testLeft: boolean,
                testRight: boolean
            ): Point | null => {
                const CONFLICT_RESOLUTION_OFFSET = 5;
                const neighboringSpaceVacant = (
                    xShift: number,
                    yShift: number
                ): boolean => {
                    for (let row = Math.floor(
                        Math.min(
                            Math.max(0, startRow + yShift),
                            this._annotationSpaceAvailability.length
                        )
                    );
                        row < Math.ceil(
                            Math.min(
                                Math.max(0, stopRow + yShift),
                                this._annotationSpaceAvailability.length
                            )
                        );
                        row++) {
                        for (let col = Math.floor(
                            Math.min(
                                Math.max(0, startCol + xShift),
                                this._annotationSpaceAvailability[0].length
                            )
                        );
                            col < Math.ceil(
                                Math.min(
                                    Math.max(0, stopCol + xShift),
                                    this._annotationSpaceAvailability[0].length
                                )
                            );
                            col++) {
                            if (!this._annotationSpaceAvailability[row][col]) {
                                return false;
                            }
                        }
                    }

                    return true;
                };

                let resolveMovingUp = testUp;
                if (testUp) {
                    resolveMovingUp = neighboringSpaceVacant(
                        0,
                        -conflict.bounds.height - CONFLICT_RESOLUTION_OFFSET
                    );

                    if (resolveMovingUp) {
                        return new Point(
                            0,
                            -conflict.bounds.height - CONFLICT_RESOLUTION_OFFSET
                        );
                    }
                }

                let resolveMovingDown = testDown;
                if (testDown) {
                    resolveMovingDown = neighboringSpaceVacant(
                        0,
                        conflict.bounds.height + CONFLICT_RESOLUTION_OFFSET
                    );

                    if (resolveMovingDown) {
                        return new Point(
                            0,
                            conflict.bounds.height + CONFLICT_RESOLUTION_OFFSET
                        );
                    }
                }

                let resolveMovingLeft = testLeft;
                if (testLeft) {
                    resolveMovingLeft = neighboringSpaceVacant(
                        -conflict.bounds.width - CONFLICT_RESOLUTION_OFFSET,
                        0
                    );

                    if (resolveMovingLeft) {
                        return new Point(
                            -conflict.bounds.width - CONFLICT_RESOLUTION_OFFSET,
                            0
                        );
                    }
                }

                let resolveMovingRight = testRight;
                if (testRight) {
                    resolveMovingRight = neighboringSpaceVacant(
                        conflict.bounds.width + CONFLICT_RESOLUTION_OFFSET,
                        0
                    );

                    if (resolveMovingRight) {
                        return new Point(
                            conflict.bounds.width + CONFLICT_RESOLUTION_OFFSET,
                            0
                        );
                    }
                }

                let resolveMovingUpLeft = testUp && testLeft;
                if (testUp && testLeft) {
                    resolveMovingUpLeft = neighboringSpaceVacant(
                        -conflict.bounds.width - CONFLICT_RESOLUTION_OFFSET,
                        -conflict.bounds.height - CONFLICT_RESOLUTION_OFFSET
                    );

                    if (resolveMovingUpLeft) {
                        return new Point(
                            -conflict.bounds.width - CONFLICT_RESOLUTION_OFFSET,
                            -conflict.bounds.height - CONFLICT_RESOLUTION_OFFSET
                        );
                    }
                }

                let resolveMovingUpRight = testUp && testRight;
                if (testUp && testRight) {
                    resolveMovingUpRight = neighboringSpaceVacant(
                        conflict.bounds.width + CONFLICT_RESOLUTION_OFFSET,
                        -conflict.bounds.height - CONFLICT_RESOLUTION_OFFSET
                    );

                    if (resolveMovingUpRight) {
                        return new Point(
                            conflict.bounds.width + CONFLICT_RESOLUTION_OFFSET,
                            -conflict.bounds.height - CONFLICT_RESOLUTION_OFFSET
                        );
                    }
                }

                let resolveMovingDownLeft = testDown && testLeft;
                if (testDown && testLeft) {
                    resolveMovingDownLeft = neighboringSpaceVacant(
                        -conflict.bounds.width - CONFLICT_RESOLUTION_OFFSET,
                        conflict.bounds.height + CONFLICT_RESOLUTION_OFFSET
                    );

                    if (resolveMovingDownLeft) {
                        return new Point(
                            -conflict.bounds.width - CONFLICT_RESOLUTION_OFFSET,
                            conflict.bounds.height + CONFLICT_RESOLUTION_OFFSET
                        );
                    }
                }

                let resolveMovingDownRight = testDown && testRight;
                if (testDown && testRight) {
                    resolveMovingDownRight = neighboringSpaceVacant(
                        conflict.bounds.width + CONFLICT_RESOLUTION_OFFSET,
                        conflict.bounds.height + CONFLICT_RESOLUTION_OFFSET
                    );

                    if (resolveMovingDownRight) {
                        return new Point(
                            conflict.bounds.width + CONFLICT_RESOLUTION_OFFSET,
                            conflict.bounds.height + CONFLICT_RESOLUTION_OFFSET
                        );
                    }
                }

                return null;
            };

            let topLeftConflict: AnnotationBaseConflict | null = null;
            let topCenterConflict: AnnotationBaseConflict | null = null;
            let topRightConflict: AnnotationBaseConflict | null = null;
            let leftCenterConflict: AnnotationBaseConflict | null = null;
            let rightCenterConflict: AnnotationBaseConflict | null = null;
            let bottomLeftConflict: AnnotationBaseConflict | null = null;
            let bottomCenterConflict: AnnotationBaseConflict | null = null;
            let bottomRightConflict: AnnotationBaseConflict | null = null;

            // Get base layer bounds relative to Pose2D container
            const baseLayerBounds = DisplayUtil.getBoundsRelative(this._baseLayer, this.container);

            // Determine anchor co-ordinates relative to base layer
            // we subtract (width / 2) from X to compute value relative to left edge versus midpoint
            // we subtract (height / 2) from Y to compute value relative to top edge versus midpoint
            const anchorBaseX: number = (anchorPoint.x + anchorOffsetX) - baseLayerBounds.x - anchorDisplay.width / 2;
            const anchorBaseY: number = (anchorPoint.y + anchorOffsetY) - baseLayerBounds.y;

            // Compute extents used to iterate across places:
            // - Center Row = Left Center and Right Center
            // - Top Row = Top Left and Top Right
            // - Bottom Row = Bottom Left and Bottom Right
            // - Center Column = Top and Bottom
            // - Left Column = Top Left, Bottom Left, Left Center
            // - Right Column = Top Rigth, Bottom right, Right Center
            const startCenterRow = anchorBaseY - annotationCard.display.height / 2;
            const stopCenterRow = anchorBaseY + annotationCard.display.height / 2;
            const startCenterColumn = anchorBaseX - annotationCard.display.width / 2 + anchorDisplay.width / 2;
            const stopCenterColumn = anchorBaseX + annotationCard.display.width / 2 + anchorDisplay.width / 2;
            const startTopRow = anchorBaseY - anchorDisplay.height / 2 - (yOffset + annotationCard.display.height);
            const stopTopRow = anchorBaseY - anchorDisplay.height / 2 - yOffset;
            const startLeftColumn = anchorBaseX - (xOffset + annotationCard.display.width);
            const stopLeftColumn = anchorBaseX - xOffset;
            const startBottomRow = anchorBaseY + anchorDisplay.height / 2 + yOffset;
            const stopBottomRow = anchorBaseY + anchorDisplay.height / 2 + yOffset + annotationCard.display.height;
            const startRightColumn = anchorBaseX + anchorDisplay.width + xOffset;
            const stopRightColumn = anchorBaseX + anchorDisplay.width + xOffset + annotationCard.display.width;

            // Top Left
            topLeftConflict = inspectRegion(
                startTopRow,
                stopTopRow,
                startLeftColumn,
                stopLeftColumn
            );

            // determine if any conflict is resolvable
            if (topLeftConflict) {
                if (DEBUG_ANNOTATION_PLACEMENT) {
                    const debugRect = new Graphics().lineStyle(1, BASE_CONFLICT_COLOR).drawRect(
                        baseLayerBounds.x + topLeftConflict.bounds.x,
                        baseLayerBounds.y + topLeftConflict.bounds.y,
                        topLeftConflict.bounds.width,
                        topLeftConflict.bounds.height
                    );
                    if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                        this._baseLayer.addChild(debugRect);
                    }
                }

                const correction = testForConflictCorrection(
                    topLeftConflict,
                    startTopRow,
                    stopTopRow,
                    startLeftColumn,
                    stopLeftColumn,
                    true,
                    false,
                    true,
                    false
                );

                if (correction) {
                    topLeftConflict.resolvable = true;
                    topLeftConflict.correction = correction;
                }
            } else if (DEBUG_ANNOTATION_PLACEMENT) {
                const debugRect = new Graphics().lineStyle(1, POSSIBLE_PLACEMENT_AVAILABILITY_COLOR).drawRect(
                    startLeftColumn + baseLayerBounds.x,
                    startTopRow + baseLayerBounds.y,
                    annotationCard.display.width,
                    annotationCard.display.height
                );
                if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                    this._baseLayer.addChild(debugRect);
                }
            }

            // Top Center
            topCenterConflict = inspectRegion(
                startTopRow,
                stopTopRow,
                startCenterColumn,
                stopCenterColumn
            );

            // determine if any conflict is resolvable
            if (topCenterConflict) {
                if (DEBUG_ANNOTATION_PLACEMENT) {
                    const debugRect = new Graphics().lineStyle(1, BASE_CONFLICT_COLOR).drawRect(
                        baseLayerBounds.x + topCenterConflict.bounds.x,
                        baseLayerBounds.y + topCenterConflict.bounds.y,
                        topCenterConflict.bounds.width,
                        topCenterConflict.bounds.height
                    );
                    if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                        this._baseLayer.addChild(debugRect);
                    }
                }

                const correction = testForConflictCorrection(
                    topCenterConflict,
                    startTopRow,
                    stopTopRow,
                    startCenterColumn,
                    stopCenterColumn,
                    true,
                    false,
                    false,
                    false
                );

                if (correction) {
                    topCenterConflict.resolvable = true;
                    topCenterConflict.correction = correction;
                }
            } else if (DEBUG_ANNOTATION_PLACEMENT) {
                const debugRect = new Graphics().lineStyle(1, POSSIBLE_PLACEMENT_AVAILABILITY_COLOR).drawRect(
                    startCenterColumn + baseLayerBounds.x,
                    startTopRow + baseLayerBounds.y,
                    annotationCard.display.width,
                    annotationCard.display.height
                );
                if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                    this._baseLayer.addChild(debugRect);
                }
            }

            // Top Right
            topRightConflict = inspectRegion(
                startTopRow,
                stopTopRow,
                startRightColumn,
                stopRightColumn
            );

            // determine if any conflict is resolvable
            if (topRightConflict) {
                if (DEBUG_ANNOTATION_PLACEMENT) {
                    const debugRect = new Graphics().lineStyle(1, BASE_CONFLICT_COLOR).drawRect(
                        baseLayerBounds.x + topRightConflict.bounds.x,
                        baseLayerBounds.y + topRightConflict.bounds.y,
                        topRightConflict.bounds.width,
                        topRightConflict.bounds.height
                    );
                    if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                        this._baseLayer.addChild(debugRect);
                    }
                }

                const correction = testForConflictCorrection(
                    topRightConflict,
                    startTopRow,
                    stopTopRow,
                    startRightColumn,
                    stopRightColumn,
                    true,
                    false,
                    false,
                    true
                );

                if (correction) {
                    topRightConflict.resolvable = true;
                    topRightConflict.correction = correction;
                }
            } else if (DEBUG_ANNOTATION_PLACEMENT) {
                const debugRect = new Graphics().lineStyle(1, POSSIBLE_PLACEMENT_AVAILABILITY_COLOR).drawRect(
                    startRightColumn + baseLayerBounds.x,
                    startTopRow + baseLayerBounds.y,
                    annotationCard.display.width,
                    annotationCard.display.height
                );
                if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                    this._baseLayer.addChild(debugRect);
                }
            }

            // Left Center
            leftCenterConflict = inspectRegion(
                startCenterRow,
                stopCenterRow,
                startLeftColumn,
                stopLeftColumn
            );

            // determine if any conflict is resolvable
            if (leftCenterConflict) {
                if (DEBUG_ANNOTATION_PLACEMENT) {
                    const debugRect = new Graphics().lineStyle(1, BASE_CONFLICT_COLOR).drawRect(
                        baseLayerBounds.x + leftCenterConflict.bounds.x,
                        baseLayerBounds.y + leftCenterConflict.bounds.y,
                        leftCenterConflict.bounds.width,
                        leftCenterConflict.bounds.height
                    );
                    if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                        this._baseLayer.addChild(debugRect);
                    }
                }

                const correction = testForConflictCorrection(
                    leftCenterConflict,
                    startCenterRow,
                    stopCenterRow,
                    startLeftColumn,
                    stopLeftColumn,
                    false,
                    false,
                    true,
                    false
                );

                if (correction) {
                    leftCenterConflict.resolvable = true;
                    leftCenterConflict.correction = correction;
                }
            } else if (DEBUG_ANNOTATION_PLACEMENT) {
                const debugRect = new Graphics().lineStyle(1, POSSIBLE_PLACEMENT_AVAILABILITY_COLOR).drawRect(
                    startLeftColumn + baseLayerBounds.x,
                    startCenterRow + baseLayerBounds.y,
                    annotationCard.display.width,
                    annotationCard.display.height
                );
                if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                    this._baseLayer.addChild(debugRect);
                }
            }

            // Right Center
            rightCenterConflict = inspectRegion(
                startCenterRow,
                stopCenterRow,
                startRightColumn,
                stopRightColumn
            );

            // determine if any conflict is resolvable
            if (rightCenterConflict) {
                if (DEBUG_ANNOTATION_PLACEMENT) {
                    const debugRect = new Graphics().lineStyle(1, BASE_CONFLICT_COLOR).drawRect(
                        baseLayerBounds.x + rightCenterConflict.bounds.x,
                        baseLayerBounds.y + rightCenterConflict.bounds.y,
                        rightCenterConflict.bounds.width,
                        rightCenterConflict.bounds.height
                    );
                    if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                        this._baseLayer.addChild(debugRect);
                    }
                }

                const correction = testForConflictCorrection(
                    rightCenterConflict,
                    startCenterRow,
                    stopCenterRow,
                    startRightColumn,
                    stopRightColumn,
                    false,
                    false,
                    false,
                    true
                );

                if (correction) {
                    rightCenterConflict.resolvable = true;
                    rightCenterConflict.correction = correction;
                }
            } else if (DEBUG_ANNOTATION_PLACEMENT) {
                const debugRect = new Graphics().lineStyle(1, POSSIBLE_PLACEMENT_AVAILABILITY_COLOR).drawRect(
                    startRightColumn + baseLayerBounds.x,
                    startCenterRow + baseLayerBounds.y,
                    annotationCard.display.width,
                    annotationCard.display.height
                );
                if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                    this._baseLayer.addChild(debugRect);
                }
            }

            // Bottom Right
            bottomRightConflict = inspectRegion(
                startBottomRow,
                stopBottomRow,
                startRightColumn,
                stopRightColumn
            );

            // determine if any conflict is resolvable
            if (bottomRightConflict) {
                if (DEBUG_ANNOTATION_PLACEMENT) {
                    const debugRect = new Graphics().lineStyle(1, BASE_CONFLICT_COLOR).drawRect(
                        baseLayerBounds.x + bottomRightConflict.bounds.x,
                        baseLayerBounds.y + bottomRightConflict.bounds.y,
                        bottomRightConflict.bounds.width,
                        bottomRightConflict.bounds.height
                    );
                    if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                        this._baseLayer.addChild(debugRect);
                    }
                }

                const correction = testForConflictCorrection(
                    bottomRightConflict,
                    startBottomRow,
                    stopBottomRow,
                    startRightColumn,
                    stopRightColumn,
                    false,
                    true,
                    false,
                    true
                );

                if (correction) {
                    bottomRightConflict.resolvable = true;
                    bottomRightConflict.correction = correction;
                }
            } else if (DEBUG_ANNOTATION_PLACEMENT) {
                const debugRect = new Graphics().lineStyle(1, POSSIBLE_PLACEMENT_AVAILABILITY_COLOR).drawRect(
                    startRightColumn + baseLayerBounds.x,
                    startBottomRow + baseLayerBounds.y,
                    annotationCard.display.width,
                    annotationCard.display.height
                );
                if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                    this._baseLayer.addChild(debugRect);
                }
            }

            // Bottom Center
            bottomCenterConflict = inspectRegion(
                startBottomRow,
                stopBottomRow,
                startCenterColumn,
                stopCenterColumn
            );

            // determine if any conflict is resolvable
            if (bottomCenterConflict) {
                if (DEBUG_ANNOTATION_PLACEMENT) {
                    const debugRect = new Graphics().lineStyle(1, BASE_CONFLICT_COLOR).drawRect(
                        baseLayerBounds.x + bottomCenterConflict.bounds.x,
                        baseLayerBounds.y + bottomCenterConflict.bounds.y,
                        bottomCenterConflict.bounds.width,
                        bottomCenterConflict.bounds.height
                    );
                    if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                        this._baseLayer.addChild(debugRect);
                    }
                }

                const correction = testForConflictCorrection(
                    bottomCenterConflict,
                    startBottomRow,
                    stopBottomRow,
                    startCenterColumn,
                    stopCenterColumn,
                    false,
                    true,
                    false,
                    false
                );

                if (correction) {
                    bottomCenterConflict.resolvable = true;
                    bottomCenterConflict.correction = correction;
                }
            } else if (DEBUG_ANNOTATION_PLACEMENT) {
                const debugRect = new Graphics().lineStyle(1, POSSIBLE_PLACEMENT_AVAILABILITY_COLOR).drawRect(
                    startCenterColumn + baseLayerBounds.x,
                    startBottomRow + baseLayerBounds.y,
                    annotationCard.display.width,
                    annotationCard.display.height
                );
                if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                    this._baseLayer.addChild(debugRect);
                }
            }

            // Bottom Left
            bottomLeftConflict = inspectRegion(
                startBottomRow,
                stopBottomRow,
                startLeftColumn,
                stopLeftColumn
            );

            // determine if any conflict is resolvable
            if (bottomLeftConflict) {
                if (DEBUG_ANNOTATION_PLACEMENT) {
                    const debugRect = new Graphics().lineStyle(1, BASE_CONFLICT_COLOR).drawRect(
                        baseLayerBounds.x + bottomLeftConflict.bounds.x,
                        baseLayerBounds.y + bottomLeftConflict.bounds.y,
                        bottomLeftConflict.bounds.width,
                        bottomLeftConflict.bounds.height
                    );
                    if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                        this._baseLayer.addChild(debugRect);
                    }
                }

                const correction = testForConflictCorrection(
                    bottomLeftConflict,
                    startBottomRow,
                    stopBottomRow,
                    startLeftColumn,
                    stopLeftColumn,
                    false,
                    true,
                    true,
                    false
                );

                if (correction) {
                    bottomLeftConflict.resolvable = true;
                    bottomLeftConflict.correction = correction;
                }
            } else if (DEBUG_ANNOTATION_PLACEMENT) {
                const debugRect = new Graphics().lineStyle(1, POSSIBLE_PLACEMENT_AVAILABILITY_COLOR).drawRect(
                    startLeftColumn + baseLayerBounds.x,
                    startBottomRow + baseLayerBounds.y,
                    annotationCard.display.width,
                    annotationCard.display.height
                );
                if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                    this._baseLayer.addChild(debugRect);
                }
            }

            return {
                [AnnotationPlacement.TOP_LEFT]: topLeftConflict,
                [AnnotationPlacement.TOP_CENTER]: topCenterConflict,
                [AnnotationPlacement.TOP_RIGHT]: topRightConflict,
                [AnnotationPlacement.LEFT_CENTER]: leftCenterConflict,
                [AnnotationPlacement.RIGHT_CENTER]: rightCenterConflict,
                [AnnotationPlacement.BOTTOM_LEFT]: bottomLeftConflict,
                [AnnotationPlacement.BOTTOM_CENTER]: bottomCenterConflict,
                [AnnotationPlacement.BOTTOM_RIGHT]: bottomRightConflict
            };
        };

        /**
         * Helper function that searches for a suitable region to attempt to place annotation
         *
         *                      top
         *     top-left   --------------- top-right
         *               |       |       |
         *               |       |       |
         *               |       |       |
         *   left-center  ---- Anchor ---  right-center
         *               |       |       |
         *               |       |       |
         *               |       |       |
         *   bottom-left  ---------------  bottom-right
         *                    bottom
         *
         * @param baseConflicts an object with the conflict (or lack thereof) at each placement region
         * @param includeCentersPosition whether to include central locations as suitable regions
         * @return relative position of existing annotion (if one exists) or null (if position is vacant)
         */
        const findProposedPosition = (
            baseConflicts: AnnotationBaseConflicts,
            includeTopLeft: boolean,
            includeTopCenter: boolean,
            includeTopRight: boolean,
            includeLeftCenter: boolean,
            includeRightCenter: boolean,
            includeBottomLeft: boolean,
            includeBottomCenter: boolean,
            includeBottomRight: boolean
        ): AnnotationPosition | null => {
            // IMPORTANT: The conditional statements that follow
            // have been ordered intentionally to
            // attempt to present the annotations in
            // this specific priority list

            if (
                includeRightCenter
                && !baseConflicts[AnnotationPlacement.RIGHT_CENTER]
            ) {
                // Place at right-center
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        anchorDisplay.width / 2 + xOffset + anchorOffsetX,
                        -annotationCard.display.height / 2 + anchorOffsetY
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.RIGHT_CENTER,
                    custom: false
                };
            }

            if (
                includeLeftCenter
                && !baseConflicts[AnnotationPlacement.LEFT_CENTER]
            ) {
                // Place at left-center
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        -anchorDisplay.width / 2 - annotationCard.display.width - xOffset + anchorOffsetX,
                        -annotationCard.display.height / 2 + anchorOffsetY
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.LEFT_CENTER,
                    custom: false
                };
            }

            if (
                includeTopCenter
                && !baseConflicts[AnnotationPlacement.TOP_CENTER]
            ) {
                // Place at top-center
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        -annotationCard.display.width / 2 + anchorDisplay.width / 2 + anchorOffsetX,
                        -anchorDisplay.height / 2 - annotationCard.display.height - yOffset + anchorOffsetY
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.TOP_CENTER,
                    custom: false
                };
            }

            if (
                includeBottomCenter
                && !baseConflicts[AnnotationPlacement.BOTTOM_CENTER]
            ) {
                // Place at bottom-center
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        -annotationCard.display.width / 2 + anchorDisplay.width / 2 + anchorOffsetX,
                        anchorDisplay.height / 2 + yOffset + anchorOffsetY
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.BOTTOM_CENTER,
                    custom: false
                };
            }

            if (
                includeTopLeft
                && !baseConflicts[AnnotationPlacement.TOP_LEFT]
            ) {
                // Place in top-left
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        -anchorDisplay.width / 2 - annotationCard.display.width - xOffset + anchorOffsetX,
                        -anchorDisplay.height / 2 - annotationCard.display.height - yOffset + anchorOffsetY
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.TOP_LEFT,
                    custom: false
                };
            }

            if (
                includeTopRight
                && !baseConflicts[AnnotationPlacement.TOP_RIGHT]
            ) {
                // Place in top-right
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        anchorDisplay.width / 2 + xOffset + anchorOffsetX,
                        -anchorDisplay.height / 2 - annotationCard.display.height - yOffset + anchorOffsetY
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.TOP_RIGHT,
                    custom: false
                };
            }

            if (
                includeBottomLeft
                && !baseConflicts[AnnotationPlacement.BOTTOM_LEFT]
            ) {
                // Place in bottom-left
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        -anchorDisplay.width / 2 - annotationCard.display.width - xOffset + anchorOffsetX,
                        anchorDisplay.height / 2 + yOffset + anchorOffsetY
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.BOTTOM_LEFT,
                    custom: false
                };
            }

            if (
                includeBottomRight
                && !baseConflicts[AnnotationPlacement.BOTTOM_RIGHT]
            ) {
                // Place in bottom-right
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        anchorDisplay.width / 2 + xOffset + anchorOffsetX,
                        anchorDisplay.height / 2 + yOffset + anchorOffsetY
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.BOTTOM_RIGHT,
                    custom: false
                };
            }

            // Shifted

            const rightCenterConflict = baseConflicts[AnnotationPlacement.RIGHT_CENTER];
            const rightCenterCorrection = rightCenterConflict?.correction;
            if (
                includeRightCenter
                && rightCenterConflict
                && rightCenterConflict?.resolvable
                && rightCenterCorrection
            ) {
                // Place in shifted right-center

                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        anchorDisplay.width / 2 + xOffset + anchorOffsetX
                        + rightCenterCorrection.x,
                        -annotationCard.display.height / 2 + anchorOffsetY
                        + rightCenterCorrection.y
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.RIGHT_CENTER,
                    custom: false
                };
            }

            const leftCenterConflict = baseConflicts[AnnotationPlacement.LEFT_CENTER];
            const leftCenterCorrection = leftCenterConflict?.correction;
            if (
                includeLeftCenter
                && leftCenterConflict
                && leftCenterConflict?.resolvable
                && leftCenterCorrection
            ) {
                // Place in shifted left-center
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        -anchorDisplay.width / 2 - annotationCard.display.width - xOffset + anchorOffsetX
                        + leftCenterCorrection.x,
                        -annotationCard.display.height / 2 + anchorOffsetY
                        + leftCenterCorrection.y
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.LEFT_CENTER,
                    custom: false
                };
            }

            const topCenterConflict = baseConflicts[AnnotationPlacement.TOP_CENTER];
            const topCenterCorrection = topCenterConflict?.correction;
            if (
                includeTopCenter
                && topCenterConflict
                && topCenterConflict?.resolvable
                && topCenterCorrection
            ) {
                // Place in shifted top-center
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        -annotationCard.display.width / 2 + anchorDisplay.width / 2 + anchorOffsetX
                        + topCenterCorrection.x,
                        -anchorDisplay.height / 2 - annotationCard.display.height - yOffset + anchorOffsetY
                        + topCenterCorrection.y
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.TOP_CENTER,
                    custom: false
                };
            }

            const bottomCenterConflict = baseConflicts[AnnotationPlacement.BOTTOM_CENTER];
            const bottomCenterCorrection = topCenterConflict?.correction;
            if (
                includeBottomCenter
                && bottomCenterConflict
                && bottomCenterConflict?.resolvable
                && bottomCenterCorrection
            ) {
                // Place in shifted bottom-center
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        -annotationCard.display.width / 2 + anchorDisplay.width / 2 + anchorOffsetX
                        + bottomCenterCorrection.x,
                        anchorDisplay.height / 2 + yOffset + anchorOffsetY
                        + bottomCenterCorrection.y
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.BOTTOM_CENTER,
                    custom: false
                };
            }

            const topLeftConflict = baseConflicts[AnnotationPlacement.TOP_LEFT];
            const topLeftCorrection = topCenterConflict?.correction;
            if (
                includeTopLeft
                && topLeftConflict
                && topLeftConflict?.resolvable
                && topLeftCorrection
            ) {
                // Place in shifted top-left
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        -anchorDisplay.width / 2 - annotationCard.display.width - xOffset
                        + anchorOffsetX + topLeftCorrection.x,
                        -anchorDisplay.height / 2 - annotationCard.display.height - yOffset
                        + anchorOffsetY + topLeftCorrection.y
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.TOP_LEFT,
                    custom: false
                };
            }

            const topRightConflict = baseConflicts[AnnotationPlacement.TOP_RIGHT];
            const topRightCorrection = topCenterConflict?.correction;
            if (
                includeTopRight
                && topRightConflict
                && topRightConflict?.resolvable
                && topRightCorrection
            ) {
                // Place in shifted top-right
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        anchorDisplay.width / 2 + xOffset + anchorOffsetX
                        + topRightCorrection.x,
                        -anchorDisplay.height / 2 - annotationCard.display.height - yOffset + anchorOffsetY
                        + topRightCorrection.y
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.TOP_RIGHT,
                    custom: false
                };
            }

            const bottomLeftConflict = baseConflicts[AnnotationPlacement.BOTTOM_LEFT];
            const bottomLeftCorrection = topCenterConflict?.correction;
            if (
                includeBottomLeft
                && bottomLeftConflict
                && bottomLeftConflict?.resolvable
                && bottomLeftCorrection
            ) {
                // Place in shifted bottom-left
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        -anchorDisplay.width / 2 - annotationCard.display.width - xOffset + anchorOffsetX
                        + bottomLeftCorrection.x,
                        anchorDisplay.height / 2 + yOffset + anchorOffsetY
                        + bottomLeftCorrection.y
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.BOTTOM_LEFT,
                    custom: false
                };
            }

            const bottomRightConflict = baseConflicts[AnnotationPlacement.BOTTOM_RIGHT];
            const bottomRightCorrection = topCenterConflict?.correction;
            if (
                includeBottomRight
                && bottomRightConflict
                && bottomRightConflict?.resolvable
                && bottomRightCorrection
            ) {
                // Place in shifted bottom-right
                return {
                    anchorIndex: currentAnchorIndex,
                    relPosition: new Point(
                        anchorDisplay.width / 2 + xOffset + anchorOffsetX
                        + bottomRightCorrection.x,
                        anchorDisplay.height / 2 + yOffset + anchorOffsetY
                        + bottomRightCorrection.y
                    ),
                    zoomLevel: this._zoomLevel,
                    placement: AnnotationPlacement.BOTTOM_RIGHT,
                    custom: false
                };
            }

            return null;
        };

        /**
         * Helper function that checks whether annotations/layers exist at a proposed position
         *
         * @param position proposed relative position (to anchor) of annotation computed from
         * the annotation's top-left corner
         * @return absolute bounds of existing annotion (if one exists) or null (if position is vacant)
         */
        const checkIfCardAtPosition = (position: AnnotationPosition | null): AnnotationPositionConflict | null => {
            if (
                position
                && position.relPosition
                && position.placement
            ) {
                const cardArray = this._zoomLevel > Pose2D.ANNOTATION_LAYER_THRESHOLD
                    ? this._layers : this._annotations;

                for (let i = 0; i < cardArray.length; i++) {
                    // Get annotation object
                    const card = cardArray[i];
                    // Annotation might have multiple positions for each range associated with it
                    for (let j = 0; j < card.positions.length; j++) {
                        const display = card.displays[j];
                        const cardRelPosition = card.positions[j].relPosition;
                        const cardAnchorPoint = new Point(
                            this._bases[card.positions[j].anchorIndex].x + this._offX,
                            this._bases[card.positions[j].anchorIndex].y + this._offY
                        );
                        const cardAbsolutePosition = new Point(
                            cardRelPosition.x + cardAnchorPoint.x,
                            cardRelPosition.y + cardAnchorPoint.y
                        );

                        const absolutePosition = new Point(
                            position.relPosition.x + anchorPoint.x,
                            position.relPosition.y + anchorPoint.y
                        );

                        // There are four cases where overlap can occur
                        if (display && ((
                            // Existing annotation behind possible position
                            // Existing annotation below possible position
                            (
                                cardAbsolutePosition.x >= absolutePosition.x
                                        && cardAbsolutePosition.x < absolutePosition.x + annotationCard.display.width
                            )
                                    && (
                                        cardAbsolutePosition.y >= absolutePosition.y
                                        && cardAbsolutePosition.y < absolutePosition.y + annotationCard.display.height
                                    )
                        )
                            || (
                                // Existing annotation behind possible position
                                // Existing annotation below possible position
                                (
                                    absolutePosition.x >= cardAbsolutePosition.x
                                    && absolutePosition.x < cardAbsolutePosition.x + display.width
                                )
                                && (
                                    cardAbsolutePosition.y >= absolutePosition.y
                                    && cardAbsolutePosition.y < absolutePosition.y + annotationCard.display.height
                                )
                            )
                            || (
                                // Existing annotation after possible position
                                // Existing annotation above possible annotation
                                (
                                    absolutePosition.x >= cardAbsolutePosition.x
                                    && absolutePosition.x < cardAbsolutePosition.x + display.width
                                )
                                && (
                                    absolutePosition.y >= cardAbsolutePosition.y
                                    && absolutePosition.y < cardAbsolutePosition.y + display.height
                                )
                            )
                            || (
                                // Existing annotation after possible position
                                // Existing annotation above possible annotation
                                (
                                    cardAbsolutePosition.x >= absolutePosition.x
                                    && cardAbsolutePosition.x < absolutePosition.x + annotationCard.display.width
                                )
                                && (
                                    absolutePosition.y >= cardAbsolutePosition.y
                                    && absolutePosition.y < cardAbsolutePosition.y + display.height
                                )
                            )
                        )) {
                            // We want absolute position so we add back anchor position
                            const positionConflict: AnnotationPositionConflict = {
                                bounds: new Rectangle(
                                    cardAbsolutePosition.x,
                                    cardAbsolutePosition.y,
                                    display.width,
                                    display.height
                                ),
                                placement: position.placement
                            };

                            if (DEBUG_ANNOTATION_PLACEMENT) {
                                const debugRect = new Graphics().lineStyle(1, ANNOTATION_CONFLICT_COLOR).drawRect(
                                    positionConflict.bounds.x,
                                    positionConflict.bounds.y,
                                    positionConflict.bounds.width,
                                    positionConflict.bounds.height
                                );
                                if (numSearchAttempts >= DEBUG_FROM_SEARCH_ATTEMPT) {
                                    this._baseLayer.addChild(debugRect);
                                }
                            }

                            return positionConflict;
                        }
                    }
                }
            }

            return null;
        };

        // Find places that are occupied by bases
        const baseConflicts = findBaseConflicts();

        // Find possible position factoring occupied places
        let proposedPosition: AnnotationPosition | null = findProposedPosition(
            baseConflicts,
            true,
            true,
            true,
            includeCenters,
            includeCenters,
            true,
            true,
            true
        );

        // Makes sure there are no annotations at proposed position
        // Will return bounds if one exists
        let annotationPositionConflict = checkIfCardAtPosition(proposedPosition);

        // Handles case if annotation exists at bounds
        // Marks region as occupied and searches next available
        // region about anchor point
        let testTopLeft = true;
        let testTopCenter = true;
        let testTopRight = true;
        let testBottomLeft = true;
        let testBottomCenter = true;
        let testBottomRight = true;
        // Accumlate all annotation overlap bounds to
        // use as anchors for recursive search
        const positionConflicts: AnnotationPositionConflict[] = [];
        while (annotationPositionConflict && proposedPosition) {
            // Store to overlap bound
            positionConflicts.push(annotationPositionConflict);

            // Update vacancy loss due to annotation occupancy
            switch (annotationPositionConflict.placement) {
                case AnnotationPlacement.TOP_LEFT:
                    testTopLeft = false;
                    break;
                case AnnotationPlacement.TOP_CENTER:
                    testTopCenter = false;
                    break;
                case AnnotationPlacement.TOP_RIGHT:
                    testTopRight = false;
                    break;
                case AnnotationPlacement.BOTTOM_LEFT:
                    testBottomLeft = false;
                    break;
                case AnnotationPlacement.BOTTOM_CENTER:
                    testBottomCenter = false;
                    break;
                case AnnotationPlacement.BOTTOM_RIGHT:
                    testBottomRight = false;
                    break;
                default:
                    break;
            }

            // Find possible position factoring updated occupied places
            proposedPosition = findProposedPosition(
                baseConflicts,
                testTopLeft,
                testTopCenter,
                testTopRight,
                false,
                false,
                testBottomLeft,
                testBottomCenter,
                testBottomRight
            );

            // Makes sure there are no annotations at new proposed position
            annotationPositionConflict = checkIfCardAtPosition(proposedPosition);
        }

        if (proposedPosition) {
            // We have an available position
            return proposedPosition.relPosition;
        } else if (positionConflicts.length > 0 && numSearchAttempts < Pose2D.ANNOTATION_PLACEMENT_ITERATION_TIMEOUT) {
            // If we still don't have a proposed position
            // Recursively search for one using each conflict annotation
            // as an anchor point
            for (const positionConflict of positionConflicts) {
                // Compute offset based on overlap placement
                const conflictOffsetX = positionConflict.bounds.x - anchorPoint.x + positionConflict.bounds.width / 2;
                const conflictOffsetY = positionConflict.bounds.y - anchorPoint.y + positionConflict.bounds.height / 2;

                const point = this.computeAnnotationPositionPoint(
                    originalAnchorIndex,
                    currentAnchorIndex,
                    anchorPoint,
                    annotationCard.display,
                    annotationCard,
                    numSearchAttempts + 1, // Increment recursive depth
                    conflictOffsetX,
                    conflictOffsetY,
                    numSearchAttempts === 0 || !includeCenters
                );

                if (point) {
                    return point;
                }
            }
        } else if (
            numSearchAttempts < Pose2D.ANNOTATION_PLACEMENT_ITERATION_TIMEOUT
            && currentAnchorIndex > 1
            && currentAnchorIndex < this._bases.length - 1
        ) {
            // We'll change the anchor index in the hopes of finding available space
            // We move in the direction with the most bases
            const increaseAnchor = this._bases.length - originalAnchorIndex > originalAnchorIndex;
            const newAnchorIndex = increaseAnchor ? currentAnchorIndex + 1 : currentAnchorIndex - 1;
            const newAnchorPoint = new Point(
                this._bases[newAnchorIndex].x + this._offX,
                this._bases[newAnchorIndex].y + this._offY
            );

            const point = this.computeAnnotationPositionPoint(
                originalAnchorIndex,
                newAnchorIndex,
                newAnchorPoint,
                annotationCard.display,
                annotationCard,
                numSearchAttempts + 1, // Increment iterative steps,
                0,
                0
            );

            if (point) {
                return point;
            }
        }

        return null;
    }

    private checkPairs(): void {
        const fullSeq = this.fullSequence;

        for (let ii = 0; ii < this._pairs.length; ii++) {
            const pi = this._pairs.pairingPartner(ii);
            if (this._pairs.isPaired(ii) && this.isPairSatisfied(ii, pi)) {
                const pairStr: number = Pose2D.getPairStrength(
                    fullSeq.nt(ii), fullSeq.nt(pi)
                );

                if (this.isAnimating) {
                    Assert.assertIsDefined(this._baseToX);
                    Assert.assertIsDefined(this._baseToY);
                    this._bases[ii].setPairing(true,
                        this._baseToX[pi] - this._baseToX[ii],
                        this._baseToY[pi] - this._baseToY[ii],
                        0.5, pairStr);
                } else {
                    this._bases[ii].setPairing(true,
                        this._bases[pi].x - this._bases[ii].x,
                        this._bases[pi].y - this._bases[ii].y,
                        0.5, pairStr);
                }
            } else {
                this._bases[ii].setPairing(false, -1, -1, 0.5, -1);
            }
        }
    }

    private updateScoreNodeVisualization(offsetChanged: boolean): void {
        if (this._scoreNodes == null) {
            this._scoreNodeHighlight.clear();
            return;
        }

        if (this.isAnimating) {
            this._scoreNodeIndex = -1;
        }

        if (this._scoreNodeIndex !== this._lastScoreNodeIndex || offsetChanged) {
            this._scoreNodeHighlight.clear();

            if (
                this._scoreNodeIndex >= 0
                && this._scoreNodes[this._scoreNodeIndex] != null
                && this._scoreNodes[this._scoreNodeIndex].baseIndices !== null
            ) {
                const origIndices = this._scoreNodes[this._scoreNodeIndex].baseIndices;
                Assert.assertIsDefined(origIndices);
                this._scoreNodeHighlight.lineStyle(0, 0, 0);
                this._scoreNodeHighlight.beginFill(0xFFFFFF, 0.22);
                const indices: number[] = origIndices.slice();

                const contour: number[] = [];
                for (let ii = 0; ii < indices.length; ii++) {
                    const p: Point = this.getBaseLoc(indices[ii]);
                    contour.push(p.x);
                    contour.push(p.y);
                }
                const triangleVerts = triangulate(contour);
                for (let ii = 0; ii < triangleVerts.length / 6; ii++) {
                    this._scoreNodeHighlight.moveTo(triangleVerts[6 * ii], triangleVerts[6 * ii + 1]);
                    this._scoreNodeHighlight.lineTo(triangleVerts[6 * ii + 2], triangleVerts[6 * ii + 3]);
                    this._scoreNodeHighlight.lineTo(triangleVerts[6 * ii + 4], triangleVerts[6 * ii + 5]);
                }
                this._scoreNodeHighlight.endFill();
            }
            this._lastScoreNodeIndex = this._scoreNodeIndex;
        }

        if (this._scoreTexts == null) return;

        for (let ii = 0; ii < this._scoreNodes.length; ii++) {
            const indices: number[] | null = this._scoreNodes[ii].baseIndices;
            Assert.assertIsDefined(indices);
            let xAvg = 0;
            let yAvg = 0;

            for (let jj = 0; jj < indices.length; jj++) {
                const p: Point = this.getBaseLoc(indices[jj]);
                xAvg += p.x;
                yAvg += p.y;
            }

            if (indices.length > 0) {
                xAvg /= indices.length;
                yAvg /= indices.length;
            }

            xAvg -= this._scoreTexts[ii].width / 2;
            yAvg -= this._scoreTexts[ii].height / 2;

            this._scoreTexts[ii].position = new Point(xAvg, yAvg);
            this._scoreTexts[ii].visible = (this._zoomLevel < 4);
            this.updateEnergyHighlight(this._scoreTexts[ii], ii, this._scoreTexts[ii].visible);
        }
    }

    public set getEnergyDelta(cb: () => number) {
        this._getEnergyDelta = cb;
    }

    private static readonly MOUSE_LOC: Point = new Point();
    private updateScoreNodeGui(): void {
        this._scoreNodeIndex = -1;

        if (this._scoreNodes != null) {
            let totalScore = 0;
            let nodeFound = false;
            let nodeLabel = '';
            let nodeScore = '';

            Assert.assertIsDefined(Flashbang.globalMouse);
            if (this._poseField.containsPoint(Flashbang.globalMouse.x, Flashbang.globalMouse.y)) {
                let mouseP: Point = new Point(0, 0);
                mouseP = mouseP.copyFrom(this.display.toLocal(Flashbang.globalMouse, undefined, Pose2D.MOUSE_LOC));
                const baseXys: Point[] = [];

                for (let ii = 0; ii < this.fullSequenceLength; ii++) {
                    baseXys.push(this.getBaseLoc(ii));
                }
                for (let ii = 0; ii < this._scoreNodes.length; ii++) {
                    const baseIndices: number[] | null = this._scoreNodes[ii].baseIndices;
                    Assert.assertIsDefined(baseIndices);
                    const nodePoints: Point[] = [];

                    for (let jj = 0; jj < baseIndices.length; jj++) {
                        nodePoints.push(baseXys[baseIndices[jj]]);
                    }

                    if (!nodeFound && Utility.isPointWithin(mouseP, nodePoints)) {
                        nodeLabel = this._scoreNodes[ii].textLabel;
                        nodeScore = this._scoreNodes[ii].textScore;
                        nodeFound = true;
                        this._scoreNodeIndex = ii;
                    }
                }
            }

            if (this._pseudoknotted && this._scoreFolder !== null) {
                totalScore = Math.round(this._scoreFolder.scoreStructures(
                    this._sequence, this._pairs.getSatisfiedPairs(this._sequence)
                ));
            } else {
                for (const scoreNode of this._scoreNodes) {
                    totalScore += scoreNode.score;
                }
            }

            let scoreLabel = 'Total';
            let scoreScore = '';
            let factor = 0;
            if ((this._molecularBindingBases != null)
                || (this._oligo != null && this._oligoMode === Pose2D.OLIGO_MODE_DIMER)
                || (this._oligos != null)) {
                const labelElems: string[] = [];
                const scoreElems: string[] = [];

                if (this._molecularBindingBases != null && this._molecularBindingBonus !== undefined) {
                    factor++;
                    if (this._moleculeIsBoundReal) {
                        labelElems.push(EnergyScoreDisplay.green('Molecule Bound'));
                        // Round to 2 decimal places
                        const bonus = Math.round(this._molecularBindingBonus * 1e2) / 1e2;
                        scoreElems.push(EnergyScoreDisplay.green(` ${bonus} kcal`));
                    } else {
                        labelElems.push(EnergyScoreDisplay.grey('Molecule Not Bound'));
                        scoreElems.push(EnergyScoreDisplay.grey(' (0 kcal)'));
                    }
                }
                if (this._oligo != null && this._oligoMode === Pose2D.OLIGO_MODE_DIMER) {
                    factor++;
                    const malus: number = this._duplexCost + Math.round(this._oligoMalus);
                    if (this._oligoPaired) {
                        labelElems.push(EnergyScoreDisplay.green('Oligo Bound'));
                        scoreElems.push(EnergyScoreDisplay.red(` ${malus.toFixed(2)} kcal`));
                    } else {
                        labelElems.push(EnergyScoreDisplay.grey('Oligo Not Bound'));
                        scoreElems.push(EnergyScoreDisplay.grey(` ${malus.toFixed(2)} kcal`));
                    }
                }
                if (this._oligos !== undefined && this._oligosOrder !== undefined) {
                    factor++;
                    if (this._oligosPaired === 0) {
                        if (this._oligos.length > 1) {
                            labelElems.push(EnergyScoreDisplay.grey('No Oligo Bound'));
                        } else {
                            labelElems.push(EnergyScoreDisplay.grey('Oligo Not Bound'));
                        }
                        scoreElems.push(EnergyScoreDisplay.grey(' (0 kcal)'));
                    } else {
                        let malus = this._duplexCost;
                        for (let ii = 0; ii < this._oligosPaired; ii++) {
                            malus += Math.round(this._oligos[this._oligosOrder[ii]].malus);
                        }
                        if (this._oligosPaired > 1) {
                            labelElems.push(EnergyScoreDisplay.green('Oligos Bound'));
                        } else {
                            labelElems.push(EnergyScoreDisplay.green('Oligo Bound'));
                        }
                        scoreElems.push(EnergyScoreDisplay.red(` ${malus.toFixed(2)} kcal`));
                    }
                }

                scoreLabel += EnergyScoreDisplay.grey(' (') + labelElems.join(', ') + EnergyScoreDisplay.grey(')');
                scoreScore = (totalScore / 100).toString() + scoreElems.join('');
            } else {
                scoreScore = `${(totalScore / 100).toString()} kcal`;
            }
            this.updateEnergyDisplaySizeLocation(factor);

            this._primaryScoreEnergyDisplay.setEnergyText(scoreLabel, scoreScore);
            this._secondaryScoreEnergyDisplay.setEnergyText(nodeLabel, nodeScore);
            this._secondaryScoreEnergyDisplay.visible = (this._showTotalEnergy && nodeFound);

            // This is because the undo stack isn't populated yet when this is run on puzzle boot/changing folders,
            // which is needed for the delta - TODO: Handle this in a less hacky way
            const attemptSetDelta = () => {
                try {
                    this._deltaScoreEnergyDisplay.setEnergyText(
                        'Natural/Target Delta',
                        `${Math.round(this._getEnergyDelta()) / 100} kcal`
                    );
                    this._deltaScoreEnergyDisplay.visible = (this._showTotalEnergy && this._scoreFolder != null);
                } catch (e) {
                    this._deltaScoreEnergyDisplay.visible = false;
                    setTimeout(attemptSetDelta, 1000);
                }
            };
            setTimeout(attemptSetDelta, 50);
        }
    }

    private updateEnergyDisplaySizeLocation(factor: number): void {
        this._primaryScoreEnergyDisplay.position = new Point(17, Pose2D.SCORES_POSITION_Y);
        this._primaryScoreEnergyDisplay.setSize(111 + factor * 59, 40);

        this._deltaScoreEnergyDisplay.position = new Point(17 + 119 + factor * 59, Pose2D.SCORES_POSITION_Y);
        this._deltaScoreEnergyDisplay.setSize(111, 40);

        this._secondaryScoreEnergyDisplay.position = new Point(17 + 119 * 2 + factor * 59, Pose2D.SCORES_POSITION_Y);
        this._secondaryScoreEnergyDisplay.setSize(111, 40);
    }

    private clearScoreTexts(): void {
        if (this._scoreTexts != null) {
            for (const scoreText of this._scoreTexts) {
                scoreText.destroy({children: true});
            }
            this._scoreTexts = null;
        }
    }

    private generateScoreNodes(): void {
        this._scoreNodes = null;
        this._scoreNodeHighlight.clear();
        this.clearEnergyHighlights();

        if (this._scoreFolder == null
            || this._sequence == null
            || this._sequence.length === 0
            || this._pairs == null
            || this._pairs.length !== this.fullSequenceLength) {
            this.clearScoreTexts();
            return;
        }

        // / JEE : It's a bit of waste to generate RNALayout twice (once here, once when drawing rna)
        // / But this is cheap, so it shouldn't matter too much
        const scoreTree: RNALayout = new RNALayout();
        scoreTree.setupTree(this.satisfiedPairs);

        const treeroot: RNATreeNode | null = scoreTree.root;
        scoreTree.scoreTree(this.fullSequence, this._scoreFolder);

        const scoreNodes: ScoreDisplayNode[] = [];
        const rootCoords: number[] = [];
        this.generateScoreNodesRecursive(treeroot, rootCoords, scoreNodes);
        this._scoreNodes = scoreNodes;

        this.clearScoreTexts();
        if (this._displayScoreTexts && !this._pseudoknotted) {
            this._scoreTexts = [];
            for (const scoreNode of this._scoreNodes) {
                const scoreText = new Sprite(BitmapManager.getTextBitmap(scoreNode.scoreString, scoreNode.scoreColor));
                scoreText.visible = false;
                this._scoreTexts.push(scoreText);
                this._energyTextLayer.addChild(scoreText);
            }
        }

        this.updateScoreNodeGui();
    }

    private generateScoreNodesRecursive(
        root: RNATreeNode | null, coords: number[] | null, nodes: ScoreDisplayNode[]
    ): void {
        if (root == null) {
            return;
        }

        if (coords != null) {
            if (root.isPair) {
                coords.push(root.indexA);
                coords.push(root.indexB);
            } else if (root.indexA >= 0) {
                coords.push(root.indexA);
                return;
            }
        }

        if (root.isPair) {
            if (root.children.length > 1) {
                throw new Error("Something's wrong with score tree");
            }

            if (root.children.length !== 0) {
                if (root.children[0].isPair) {
                    const childCoords = [];

                    childCoords.push(root.indexA);
                    childCoords.push(root.indexB);

                    childCoords.push(root.children[0].indexB);
                    childCoords.push(root.children[0].indexA);

                    const newnode = new ScoreDisplayNode();
                    nodes.push(newnode);
                    newnode.setType(ScoreDisplayNodeType.STACK, childCoords, root.score);

                    this.generateScoreNodesRecursive(root.children[0], null, nodes);
                } else {
                    const childCoords = [];

                    childCoords.push(root.indexB);
                    childCoords.push(root.indexA);

                    this.generateScoreNodesRecursive(root.children[0], childCoords, nodes);
                }
            }
        } else {
            for (const child of root.children) {
                this.generateScoreNodesRecursive(child, coords, nodes);
            }

            if (coords != null) {
                const newnode = new ScoreDisplayNode();
                nodes.push(newnode);

                newnode.setType(ScoreDisplayNodeType.LOOP, coords, root.score);
            }
        }
    }

    private isEditable(seqnum: number): boolean {
        if (this._editableIndices != null) {
            const inList: boolean = (this._editableIndices.indexOf(seqnum) !== -1);
            return this._editable ? inList : !inList;
        } else {
            return this._editable;
        }
    }

    private createBase(): Base {
        const base: Base = new Base(RNABase.GUANINE);
        this.addObject(base, this._baseLayer);
        this._bases.push(base);
        return base;
    }

    private isNucleotidePartOfSequence(index: number) {
        return index < this.fullSequence.length && this._bases[index].type !== RNABase.CUT;
    }

    private static createDefaultLocks(sequenceLength: number): boolean[] {
        const locks: boolean[] = new Array<boolean>(sequenceLength);
        for (let ii = 0; ii < sequenceLength; ++ii) {
            locks[ii] = false;
        }
        return locks;
    }

    private static getPairStrength(s1: number, s2: number): number {
        if (Pose2D.isPair(s1, s2, RNABase.ADENINE, RNABase.URACIL)) {
            return 2;
        } else if (Pose2D.isPair(s1, s2, RNABase.GUANINE, RNABase.URACIL)) {
            return 1;
        } else if (Pose2D.isPair(s1, s2, RNABase.GUANINE, RNABase.CYTOSINE)) {
            return 3;
        } else {
            return -1;
        }
    }

    private static isPair(s1: number, s2: number, type1: number, type2: number): boolean {
        return (s1 === type1 && s2 === type2) || (s1 === type2 && s2 === type1);
    }

    private readonly _baseLayer: Container = new Container();
    private readonly _poseField: PoseField;

    private _width: number = 0;
    private _height: number = 0;

    // Array of sequence/pairs
    private _sequence: Sequence = new Sequence([]);
    private _mutatedSequence: Sequence | null;
    private _pairs: SecStruct = new SecStruct();
    private _targetPairs: SecStruct = new SecStruct();
    private _pseudoknotPairs: SecStruct = new SecStruct();
    private _bases: Base[] = [];
    private _locks: boolean[] | undefined = [];
    private _forcedStruct: number[] | null = [];
    private _designStruct: boolean[] = [];
    private _bindingSite: boolean[] | null;
    private _molecularBindingBases: BaseGlow[] | null = null;
    private _molecularBindingPairs: number[] = [];
    private _molecule: Molecule | null= null;
    private _moleculeIsBound: boolean = false;
    private _moleculeIsBoundReal: boolean = false;
    private _molecularBindingBonus: number | undefined = 0;
    private _moleculeTargetPairs: SecStruct | null;
    private _shiftLimit: number;
    private _customLayout: Array<[number, number] | [null, null]> | undefined = undefined;
    private _customLayoutChanged: boolean = false;
    private _pseudoknotted: boolean = false;

    // Oligos
    private _oligo: number[] | null = null;
    private _oligoMode: number = Pose2D.OLIGO_MODE_DIMER;
    private _oligoName: string | null = null;
    private _duplexCost: number = 4.1; // total for all strands
    private _oligoMalus: number = 0; // concentration related penalty
    private _oligoBases: BaseGlow[] | null = null; // for glows
    private _oligoPaired: boolean = false;

    // Multistrands
    private _oligos: Oligo[] | undefined = undefined;
    private _oligosOrder: number[] | undefined = undefined;
    private _prevOligosOrder: number[] | undefined;
    private _oligosPaired: number = 0;
    private _strandLabel: TextBalloon;

    private _barcodes: number[];
    private _moleculeLayer: Container;
    private _energyTextLayer: Container;

    private _coloring: boolean = false;
    private _currentColor: number = RNABase.URACIL;
    private _lastColoredIndex: number;
    private _lockUpdated: boolean;
    private _bindingSiteUpdated: boolean;
    private _designStructUpdated: boolean;

    private _currentArrangementTool: Layout = Layout.MOVE;

    // Rope connecting bases for crazy user-defined layouts
    private _baseRope: BaseRope;

    // lines connecting pseudoknotted BPs
    private _pseudoknotLines: PseudoknotLines;

    // Scripted painters
    private _dynPaintColors: number[] = [];
    private _dynPaintTools: Booster[] = [];

    // Is this pose editable?
    private _editable: boolean;
    private _editableIndices: number[] | null = null;

    // Pointer to callback function to be called after change in pose
    private _poseEditCallback: (() => void) | null = null;
    private _trackMovesCallback: ((count: number, moves: Move[]) => void) | null = null;
    private _addBaseCallback: (parenthesis: string | null, op: PuzzleEditOp | null, index: number) => void;
    private _startMousedownCallback: PoseMouseDownCallback;
    private _mouseDownAltKey: boolean = false;

    // Pointer to function that needs to be called in a GameMode to have access to appropriate state
    private _getEnergyDelta: () => number;

    private _lettermode: boolean = false;
    private _displayScoreTexts: boolean;

    private _redraw: boolean = true;

    // Time which we sampled bases to animate last time;
    private lastSampledTime: number = -1;

    // Pose position offset
    private _offX: number = 0;
    private _offY: number = 0;
    private _prevOffsetX: number = 0;
    private _prevOffsetY: number = 0;
    private _offsetTranslating: boolean;
    private _startOffsetX: number;
    private _startOffsetY: number;
    private _endOffsetX: number;
    private _endOffsetY: number;

    // For base moving animation
    private _baseFromX: number[] | null;
    private _baseFromY: number[] | null;
    private _baseToX: number[] | null;
    private _baseToY: number[] | null;
    private _foldStartTime: number;
    private _foldDuration: number;
    private _paintCursor: PaintCursor;
    private _baseRotationDirectionSign: number[];

    private _zoomLevel: number = 0;
    private _desiredAngle: number = 0;

    // Is explosion animation on going?
    private _isExploding: boolean = false;
    private _explosionStartTime: number = -1;
    private _explosionRays: LightRay[];
    private _origOffsetX: number;
    private _origOffsetY: number;

    private _onExplosionComplete: (() => void) | null;

    // Selection box
    private _selectionHighlightBox: HighlightBox;
    private _restrictedHighlightBox: HighlightBox;
    private _highlightRestricted: boolean = false;
    private _unstableHighlightBox: HighlightBox;
    private _forcedHighlightBox: HighlightBox;
    private _userDefinedHighlightBox: HighlightBox;
    private _shiftHighlightBox: HighlightBox;
    private _shiftStart: number = -1;
    private _shiftEnd: number = -1;

    // For praising stacks
    private _praiseQueue: number[] = [];
    private _praiseSeq: number[] = [];

    // Score display nodes
    private _scoreNodes: ScoreDisplayNode[] | null;
    private _scoreTexts: Sprite[] | null;
    private _scoreFolder: Folder | null;
    private _scoreNodeIndex: number = -1;
    private _lastScoreNodeIndex: number = -1;
    private _scoreNodeHighlight: Graphics;

    // New Score Display panels
    private _primaryScoreEnergyDisplay: EnergyScoreDisplay;
    private _secondaryScoreEnergyDisplay: EnergyScoreDisplay;
    private _deltaScoreEnergyDisplay: EnergyScoreDisplay;
    private _showTotalEnergy: boolean = true;

    // Explosion Factor (RNALayout pairSpace multiplier)
    private _explosionFactor: number = 1;
    private _explosionFactorPanel: ExplosionFactorPanel;

    // For tracking a base
    private _cursorIndex: number | null = 0;
    private _cursorBox: Graphics | null = null;
    private _lastShiftedIndex: number = -1;
    private _lastShiftedCommand: number = -1;

    // Rendering mode
    private _numberingMode: boolean = false;
    private _showBaseRope: boolean = false;
    private _showPseudoknots: boolean = false;
    private _simpleGraphicsMods: boolean = false;
    private _annotationModeActive: boolean = false;

    // customNumbering
    private _customNumbering: (number | null)[] | undefined = undefined;

    // Last exp paint data
    private _expPainter: ExpPainter | null = null;
    private _expMid: number = 0;
    private _expHi: number = 0;
    private _expContinuous: boolean = false;
    private _expExtendedScale: boolean = false;

    private _anchoredObjects: RNAAnchorObject[] = [];
    private _highlightEnergyText: boolean = false;
    private _energyHighlights: SceneObject[] = [];

    private _showNucleotideRange: boolean = false;

    // Annotations
    private _annotations: AnnotationDisplayObject[] = [];
    private _layers: AnnotationDisplayObject[] = [];
    private _annotationGraph: AnnotationGraph;
    private _annotationSpaceAvailability: boolean[][] = [];
    private annotationCanvas: Graphics;
    private _annotationRanges: AnnotationRange[] = [];
    private _editingAnnotation: boolean = false;
    private _annotationHighlightBox: HighlightBox;
    private _puzzleAnnotationsEditable: boolean = false;
    private _resetAnnotationPositions: boolean = false;
    private _ignoreCustomAnnotationPositions: boolean = false;
    public movingAnnotation: boolean = false;

    /*
     * NEW HIGHLIGHT.
     *  - Input: List of nucleotides that we wish to highlight.
     *  - Unhighlighted Nucleotides: Draw at 65% opacity.
     *  - Highlight Nucleotides: Brighten glow around the nucleotide.
     */
    private _allNewHighlights: RNAHighlightState[] = [];

    private static readonly P: Point = new Point();

    private static readonly ANNOTATION_UNHIGHLIGHTED_OPACITY = 0.5;
    private static readonly DEFAULT_ANNOTATION_SHIFT = 15;
    private static readonly ANNOTATION_PLACEMENT_ITERATION_TIMEOUT = 20;
    private static readonly ANNOTATION_LAYER_THRESHOLD = 1;
}

export interface Oligo {
    malus: number;
    name?: string;
    sequence: number[];
}

export class RNAHighlightState {
    public nuc: number[] | null = null; // nucleotides
    public isOn: boolean = false;
}
