/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.bio.gui.glyph;

import java.awt.Shape;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.util.List;


/**
 * Utils class used by Glyph implementors. Creates Shapes from x and y coordinates
 *
 * @author Mark Southern
 * @since 1.5
 */
public class GlyphUtils {
    private GlyphUtils() {
    }

    public static Shape getShape(List x, List y) {
        GeneralPath p = new GeneralPath();

        for (int i = 0; i < Math.min(x.size(), y.size()); i++) {
            int ix = (( Integer ) x.get(0)).intValue();
            int iy = (( Integer ) y.get(0)).intValue();

            if (i == 0) {
                p.moveTo(ix, iy);
            } else {
                p.lineTo(ix, iy);
            }
        }

        p.closePath();

        return p;
    }

    public static Shape getShape(List points) {
        GeneralPath p = new GeneralPath();

        for (int i = 0; i < points.size(); i++) {
            Point2D.Float p2f = ( Point2D.Float ) points.get(i);

            if (i == 0) {
                p.moveTo(p2f.x, p2f.y);
            } else {
                p.lineTo(p2f.x, p2f.y);
            }
        }

        p.closePath();

        return p;
    }

    public static Shape getShape(int[] x, int[] y) {
        GeneralPath p = new GeneralPath();
        p.moveTo(x[ 0 ], y[ 0 ]);

        for (int i = 1; i < Math.min(x.length, y.length); i++) {
            p.lineTo(x[ i ], y[ i ]);
        }

        p.closePath();

        return p;
    }
}
