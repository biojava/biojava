/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.storage;

import java.util.List;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;

/**
 * This is a common method that can be used across multiple storage/proxy implementations to
 * handle Negative strand and other interesting elements of sequence data.
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class SequenceAsStringHelper<C extends Compound> {

    /**
     *
     * @param parsedCompounds
     * @param compoundSet
     * @param bioBegin
     * @param bioEnd
     * @param strand
     * @return
     */
    public String getSequenceAsString(List<C> parsedCompounds, CompoundSet<C> compoundSet, Integer bioBegin, Integer bioEnd, Strand strand) {
        // TODO Optimise/cache.
        if(parsedCompounds.size() == 0)
            return "";
        StringBuilder builder = new StringBuilder();
        if (strand.equals(Strand.NEGATIVE)) {
            //we expect bioBegin to be bigger but could have circular case
            if (bioBegin <= bioEnd) {
                for (int index = bioEnd - 1; index >= bioBegin - 1; index--) {
                    C compound = parsedCompounds.get(index);
                    builder.append(compoundSet.getStringForCompound(compound));
                }
            }else{
                //go to 0 and the up
                for (int index = bioBegin - 1; index >= 0; index--) {
                    C compound = parsedCompounds.get(index);
                    builder.append(compoundSet.getStringForCompound(compound));
                }
                
                for (int index = parsedCompounds.size() - 1; index >= bioEnd - 1; index--) {
                    C compound = parsedCompounds.get(index);
                    builder.append(compoundSet.getStringForCompound(compound));
                }
            }
        } else {
            if (bioBegin <= bioEnd) {
                for (int index = bioBegin - 1; index <= bioEnd - 1 ; index++) {
                    C compound = parsedCompounds.get(index);
                    builder.append(compoundSet.getStringForCompound(compound));
                }
            }else{
                //go to 0 and the up
                for (int index = bioBegin - 1; index <=  parsedCompounds.size() - 1; index++) {
                    C compound = parsedCompounds.get(index);
                    builder.append(compoundSet.getStringForCompound(compound));
                }

                for (int index = 0; index <= bioEnd - 1; index++) {
                    C compound = parsedCompounds.get(index);
                    builder.append(compoundSet.getStringForCompound(compound));
                }
            }


        }

        return builder.toString();
    }
}
