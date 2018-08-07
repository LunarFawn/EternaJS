import {DisplayObject, Graphics, Matrix} from "pixi.js";
import {Flashbang} from "../core/Flashbang";
import {GameObject} from "../core/GameObject";
import {MatrixUtil} from "../util/MatrixUtil";

/**
 * Wraps an HTML element that lives in the DOM and is drawn on top of the PIXI canvas.
 * Contains a "dummy" Container DisplayObject that mirrors the element's transform.
 */
export abstract class DOMObject<T extends HTMLElement> extends GameObject {
    protected constructor(domParentID: string, obj: T) {
        super();

        this._obj = obj;
        this._obj.style.position = "absolute";
        this._obj.style.transformOrigin = "0 0";

        // Set the initial opacity to 0 so that the object will be hidden
        // until the first postrender event. This prevents it from flickering
        // briefly on the frame it's added.
        this._obj.style.opacity = "0";

        this._domParent = document.getElementById(domParentID);
    }

    public get display(): DisplayObject {
        return this._dummyDisp;
    }

    protected added(): void {
        super.added();
        this._domParent.appendChild(this._obj);
        this.onSizeChanged();

        // Update the HTML element's transform during the PIXI postrender event -
        // this is the point where the dummy display object's transform will be up to date.
        Flashbang.pixi.renderer.addListener("postrender", this.updateElementProperties, this);
    }

    protected dispose(): void {
        this._domParent.removeChild(this._obj);
        Flashbang.pixi.renderer.removeListener("postrender", this.updateElementProperties, this);

        super.dispose();
    }

    protected updateElementProperties(): void {
        let m = this.display.worldTransform;
        if (!MatrixUtil.equals(this._lastTransform, m)) {
            this._obj.style.transform = `matrix(${m.a}, ${m.b}, ${m.c}, ${m.d}, ${m.tx}, ${m.ty})`;
            m.copy(this._lastTransform);
        }

        this._obj.style.visibility = this.display.worldVisible ? "visible" : "hidden";
        this._obj.style.opacity = this.display.worldAlpha.toString();
    }

    /**
     * Updates the dummy display object's bounds to match that of the DOM object.
     * Subclasses should call this when the DOM object's size has changed.
     */
    protected onSizeChanged(): void {
        if (this.isLiveObject) {
            let transfom: string = this._obj.style.transform;
            this._obj.style.transform = null;

            let r = this._obj.getBoundingClientRect();
            this._dummyDisp.clear()
                .beginFill(0x0, 0)
                .drawRect(0, 0, r.width, r.height)
                .endFill();

            this._obj.style.transform = transfom;
        }
    }

    /**
     * Applies the given style to the DOM object and all children who do not already have the given style property set.
     * This will not overrwrite existing properties, unless replaceIfExists is true.
     *
     * If elementNames is non-null, the style will only be applied to elements with the given names.
     */
    protected static applyStyleRecursive(element: HTMLElement, name: string, value: string,
        replaceIfExists: boolean = false, elementNames: string[] = null): void {
        let apply = true;
        if (elementNames != null) {
            apply = false;
            let thisName = element.nodeName.toUpperCase();
            for (let allowedName of elementNames) {
                if (allowedName.toUpperCase() == thisName) {
                    apply = true;
                    break;
                }
            }
        }

        if (apply && !replaceIfExists) {
            let cur = element.style.getPropertyValue(name);
            if (cur != null && cur.length > 0) {
                apply = false;
            }
        }

        if (apply) {
            element.style.setProperty(name, value);
        }

        for (let ii = 0; ii < element.children.length; ++ii) {
            let child = <HTMLElement> (element.children[ii] as any);
            if (child.accessKey !== undefined) {
                this.applyStyleRecursive(child, name, value, replaceIfExists, elementNames);
            }
        }
    }

    protected static sizeToString(size: number): string {
        return `${size}px`;
    }

    protected static stringToSize(value: string): number {
        let idx = value.indexOf("px");
        if (idx >= 0) {
            value = value.substr(0, idx);
        }
        let size: number = Number(value);
        return !isNaN(size) ? size : 0;
    }

    protected readonly _dummyDisp: Graphics = new Graphics();
    protected readonly _domParent: HTMLElement;
    protected readonly _obj: T;

    private readonly _lastTransform: Matrix = new Matrix();
}
