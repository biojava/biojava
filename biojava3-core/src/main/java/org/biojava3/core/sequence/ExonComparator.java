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
 * Created on DATE
 *
 */

package org.biojava3.core.sequence;

import java.util.Comparator;


/**
 * Sort Exon where it is a little confusing if exons shoud always be ordered left to right
 * where a negative stranded gene should go the other direction. Need to think about this?
 * @author Scooter Willis <willishf at gmail dot com>
 */
    public class ExonComparator implements Comparator<ExonSequence>{

     
        public int compare(ExonSequence o1, ExonSequence o2) {

            return o1.getBioBegin() - o2.getBioBegin();
        }

    }
