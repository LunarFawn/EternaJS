import {
    Graphics, Text
} from 'pixi.js';
import {DropShadowFilter} from '@pixi/filter-drop-shadow';
import Fonts from 'eterna/util/Fonts';
import BaseGamePanel from './BaseGamePanel';

export enum GamePanelType {
    NORMAL, INVISIBLE
}

interface GamePanelProps {
    type?: GamePanelType;
    alpha?: number;
    color?: number;
    borderAlpha?: number;
    borderColor?: number;
    dropShadow?: boolean;
    borderRadius?: number;
    borderThickness?: number;
}

export default class GamePanel extends BaseGamePanel {
    constructor(props: GamePanelProps) {
        super();

        const type = props.type || GamePanelType.NORMAL;
        const alpha = props.alpha !== undefined ? props.alpha : 0.07;
        const color = props.color ?? 0xffffff;
        const borderAlpha = props.borderAlpha !== undefined ? props.borderAlpha : 0.0;
        const borderColor = props.borderColor !== undefined ? props.borderColor : 0;
        const dropShadow = props.dropShadow || false;
        const borderRadius = props.borderRadius !== undefined ? props.borderRadius : 5;
        const borderThickness = props.borderThickness !== undefined
            ? props.borderThickness : GamePanel.DEFAULT_BORDER_THICKNESS;

        // Clicks should not pass through the panel
        this.pointerDown.connect((e) => {
            e.stopPropagation();
        });

        this._background = new Graphics();
        if (dropShadow) {
            this._background.filters = [new DropShadowFilter()];
        }
        this.container.addChild(this._background);

        this.setup(
            type,
            alpha,
            color,
            borderAlpha,
            borderColor,
            borderRadius,
            borderThickness
        );
    }

    public setup(
        type: GamePanelType,
        alpha: number,
        color: number,
        borderAlpha: number,
        borderColor: number,
        borderRadius: number = 0,
        borderThickness: number = 0
    ): void {
        this._type = type;
        this._alpha = alpha;
        this._color = color;
        this._borderAlpha = borderAlpha;
        this._borderColor = borderColor;
        this._borderRadius = borderRadius;
        this._borderThickness = borderThickness;
        this.updateView();
    }

    public setSize(width: number, height: number): void {
        this._width = width;
        this._height = height;

        this.updateView();
    }

    public setBorderThickness(thickness: number = GamePanel.DEFAULT_BORDER_THICKNESS): void {
        this._borderThickness = thickness;
        this.updateView();
    }

    public set color(color: number) {
        this._color = color;
        this.updateView();
    }

    public set alpha(alpha: number) {
        this._alpha = alpha;
        this.updateView();
    }

    public set title(title: string) {
        this._title = title;
        this.updateView();
    }

    public get titleHeight(): number {
        return this._title == null ? 0 : 35;
    }

    public get width(): number {
        return this._width;
    }

    public get height(): number {
        return this._height;
    }

    protected updateView(): void {
        this._background?.clear();

        if (this._width <= 0 || this._height <= 0) {
            return;
        }

        if (this._type === GamePanelType.INVISIBLE) {
            this._background.beginFill(0x000000, 0);
            this._background.drawRect(0, 0, this._width, this._height);
            this._background.endFill();
        } else {
            this._background.lineStyle(this._borderThickness, this._borderColor, this._borderAlpha);
            this._background.beginFill(this._color, this._alpha);
            this._background.drawRoundedRect(0, 0, this._width, this._height, this._borderRadius);
            this._background.endFill();

            if (this._title !== null) {
                if (this._titleText == null) {
                    this._titleText = Fonts.std().bold().fontSize(16).color(0xffffff)
                        .build();
                    this.container.addChild(this._titleText);
                }

                this._titleText.text = this._title.toUpperCase();
                this._titleText.position.set(
                    (this._width - this._titleText.width) * 0.5,
                    (this.titleHeight - this._titleText.height) * 0.5
                );
                this._background.beginFill(this._borderColor, this._borderAlpha);
                this._background
                    .moveTo(0, 35)
                    .lineTo(0, 5)
                    .arcTo(0, 0, 5, 0, 5)
                    .lineTo(this._width - 5, 0)
                    .arcTo(this._width, 0, this._width, 5, 5)
                    .lineTo(this._width, 35);
                this._background.endFill();
            }
        }
    }

    protected readonly _background: Graphics;

    protected _type: GamePanelType;

    protected _alpha: number = 0;
    protected _color: number = 0;
    protected _borderAlpha: number = 0;
    protected _borderColor: number = 0;
    protected _borderRadius: number = 5;
    protected _borderThickness: number = 0;
    protected _title: string | null = null;
    protected _titleText: Text | null = null;

    protected _width: number = 0;
    protected _height: number = 0;

    private static DEFAULT_BORDER_THICKNESS: number = 1.5;
}
