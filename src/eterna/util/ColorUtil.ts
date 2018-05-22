export class ColorUtil {
    /**
     * Blends two colors to create a new one.
     * @param c1 the first color
     * @param c2 the second color
     * @param blendFactor the percent contribution that c1 makes to the blended color (a value between [0, 1]).
     * (c2's contribution will be (1-blendFactor)).
     */
    public static blend (c1: number, c2: number, blendFactor: number): number {
        const r1: number = (c1 >> 16) & 0xff;
        const g1: number = (c1 >> 8) & 0xff;
        const b1: number = (c1 & 0xff);

        const r2: number = (c2 >> 16) & 0xff;
        const g2: number = (c2 >> 8) & 0xff;
        const b2: number = (c2 & 0xff);

        const blendInv: number = 1 - blendFactor;

        const rOut: number = (r1 * blendFactor) + (r2 * blendInv);
        const gOut: number = (g1 * blendFactor) + (g2 * blendInv);
        const bOut: number = (b1 * blendFactor) + (b2 * blendInv);

        return ((rOut & 0xff) << 16) | ((gOut & 0xff) << 8) | (bOut & 0xff);
    }

    /** Returns the 8-bit red component of a 24-bit color. The value will be in [0,255] */
    public static getRed (color :number) :number {
        return (color >> 16) & 0xff;
    }

    /** Returns the 8-bit green component of a 24-bit color. */
    public static getGreen (color :number) :number {
        return (color >> 8) & 0xff;
    }

    /** Returns the 8-bit blue component of a 24-bit color. */
    public static getBlue (color :number) :number {
        return color & 0xff;
    }
}