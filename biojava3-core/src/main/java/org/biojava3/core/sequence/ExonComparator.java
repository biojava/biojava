/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence;

import java.util.Comparator;


/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
    public class ExonComparator implements Comparator<ExonSequence>{

        @Override
        public int compare(ExonSequence o1, ExonSequence o2) {

            return o1.getBioBegin() - o2.getBioBegin();
        }

    }
