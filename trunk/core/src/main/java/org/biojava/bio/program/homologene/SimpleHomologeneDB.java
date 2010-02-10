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

package org.biojava.bio.program.homologene;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.biojava.utils.ChangeVetoException;

public class SimpleHomologeneDB 
    extends SimpleOrthoPairCollection 
    implements HomologeneDB
{

    // every Orthologue is stored in an OrthologueSet delegate
    private OrthologueSet orthologues = new SimpleOrthologueSet();

    // orthologies are also stored in a set
    private Set orthologySet = new HashSet();

    // indices
    private Map orthologyByTaxonID = new HashMap();
    private Map orthologyBySimilarityType = new HashMap();

    public Orthologue createOrthologue(Taxon taxon, String locusID, String homologeneID, String accession)
        throws ChangeVetoException
    {
        // create the Orthologue
        Orthologue newOrthologue = new SimpleOrthologue(taxon, locusID, homologeneID, accession);

        orthologues.addOrthologue(newOrthologue);
        return newOrthologue;
    }

    public Orthologue createOrthologue(int taxonID, String locusID, String homologeneID, String accession)
        throws IllegalArgumentException, ChangeVetoException
    {
        // create the Orthologue
        Orthologue newOrthologue = new SimpleOrthologue(taxonID, locusID, homologeneID, accession);

        orthologues.addOrthologue(newOrthologue);
        return newOrthologue;
   }

    public Orthologue getOrthologue(String homologeneID)
    {
        return orthologues.getOrthologue(homologeneID);
    }

    public OrthoPair createOrthoPair(Orthologue first, Orthologue second, SimilarityType type, double percentIdentity)
    {
        OrthoPair newOrthoPair = new SimpleOrthoPair(first, second, type, percentIdentity);

        // index it
        indexByTaxonID(first.getTaxonID(), newOrthoPair);
        indexByTaxonID(second.getTaxonID(), newOrthoPair);
        indexBySimilarityType(type, newOrthoPair);

        orthologySet.add(newOrthoPair);

        return newOrthoPair;
    }

    // should implement a uniqueness check here later!!!!

    public OrthoPair createOrthoPair(Orthologue first, Orthologue second, String ref)
    {
        OrthoPair newOrthoPair = new SimpleOrthoPair(first, second, ref);

        // index it
        indexByTaxonID(first.getTaxonID(), newOrthoPair);
        indexByTaxonID(second.getTaxonID(), newOrthoPair);
        indexBySimilarityType(SimilarityType.CURATED, newOrthoPair);

        orthologySet.add(newOrthoPair);

        return newOrthoPair;
    }

    public OrthoPairSet createOrthoPairSet()
    {
        OrthoPairSet newGroup = new SimpleOrthoPairSet();
        groups.add(newGroup);        

        return newGroup;
    }

    public OrthoPairCollection getOrthoPairSets()
    {
        return new SimpleOrthoPairCollection(Collections.unmodifiableSet(groups));
    }

    private void indexByTaxonID(int taxonID, OrthoPair orthology)
    {
        Integer taxonIDIndex = new Integer(taxonID);
        Set indexSet = (Set) orthologyByTaxonID.get(taxonIDIndex);

        if (indexSet == null) {
            indexSet = new HashSet();
            orthologyByTaxonID.put(taxonIDIndex, indexSet);
        }

        indexSet.add(orthology);
    }

    private void indexBySimilarityType(SimilarityType type, OrthoPair orthology)
    {
        Set indexSet = (Set) orthologyBySimilarityType.get(type);

        if (indexSet == null) {
            indexSet = new HashSet();
            orthologyByTaxonID.put(type, indexSet);
        }

        indexSet.add(orthology);
    }
}
