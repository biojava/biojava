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

import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;


/**
 * The Glyph interface for painting a shape within bounds
 *
 * @author Mark Southern
 * @since 1.5
 */
public interface Glyph {
    void setBounds(Rectangle2D.Float r);

    Rectangle2D.Float getBounds();

    void render(Graphics2D g);
}
