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

package org.biojava.bio.program.gff;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Comparator;

class GFFComparator implements Comparator {
    public int compare(Object o1, Object o2) {
        GFFRecord r1 = (GFFRecord) o1;
        GFFRecord r2 = (GFFRecord) o2;
        
        {
            int dif = r1.getSeqName().compareTo(r2.getSeqName());
            if (dif != 0) {
                return dif;
            }
        }
        {
            int dif = r1.getStart() - r2.getStart();
            if (dif != 0) {
                return dif;
            }
        }
        {
            int dif = r1.getEnd() - r2.getEnd();
            if (dif != 0) {
                return dif;
            }
        }
        {
            int dif = r1.getFeature().compareTo(r2.getFeature());
            if (dif != 0) {
                return dif;
            }
        }
        {
            int dif = r1.getSource().compareTo(r2.getSource());
            if (dif != 0) {
                return dif;
            }
        }
        return stringify(r1).compareTo(stringify(r2));
    }
    
    private String stringify(GFFRecord r) {
        StringWriter sw = new StringWriter();
        GFFWriter gffw = new GFFWriter(new PrintWriter(sw));
        gffw.recordLine(r);
        gffw.endDocument();
        return sw.toString();
    }
}
            
