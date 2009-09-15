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
package org.biojava.bio.gui.sequence;

import org.biojava.bio.gui.glyph.ArrowGlyph;
import org.biojava.bio.gui.glyph.HelixGlyph;
import org.biojava.bio.gui.glyph.TurnGlyph;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.utils.ChangeVetoException;


/**
 * A GlyphRenderer subclass that specificatlly handles Features pertaining to Secondary Structure
 * (Helices, Turns and Strands).
 *
 * @author Mark Southern
 * @see GlyphFeatureRenderer
 * @since 1.5
 */
public class SecondaryStructureFeatureRenderer extends GlyphFeatureRenderer {
    public SecondaryStructureFeatureRenderer() throws ChangeVetoException {
        FeatureFilter ff = new FeatureFilter.ByType("HELIX");
        addFilterAndGlyph(ff, new HelixGlyph());
        ff = new FeatureFilter.ByType("STRAND");
        addFilterAndGlyph(ff, new ArrowGlyph());
        ff = new FeatureFilter.ByType("TURN");
        addFilterAndGlyph(ff, new TurnGlyph());
    }
}
