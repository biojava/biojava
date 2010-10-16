/**
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


public interface OrthologueFilter
{
    public boolean accept(Orthologue ortho);

    public final static class AcceptAll implements OrthologueFilter
    {
        public boolean accept(Orthologue ortho) { return true; }
    }

    public final static class Not implements OrthologueFilter
    {
        OrthologueFilter a;

        public Not(OrthologueFilter a)
        {
            this.a = a;
        }

        public boolean accept(Orthologue ortho)
        {
            return !a.accept(ortho);
        }
    }

    public final static class Or implements OrthologueFilter
    {
        OrthologueFilter a;
        OrthologueFilter b;

        public Or(OrthologueFilter a, OrthologueFilter b)
        {
            this.a = a;
            this.b = b;
        }

        public boolean accept(Orthologue ortho)
        {
            return ((a.accept(ortho)) || (b.accept(ortho)));
        }
    }

    public final static class And implements OrthologueFilter
    {
        OrthologueFilter a;
        OrthologueFilter b;

        public And(OrthologueFilter a, OrthologueFilter b)
        {
            this.a = a;
            this.b = b;
        }

        public boolean accept(Orthologue ortho)
        {
            return ((a.accept(ortho)) && (b.accept(ortho)));
        }
    }

    public final static class Xor implements OrthologueFilter
    {
        OrthologueFilter a;
        OrthologueFilter b;

        public Xor(OrthologueFilter a, OrthologueFilter b)
        {
            this.a = a;
            this.b = b;
        }

        public boolean accept(Orthologue ortho)
        {
            return ((a.accept(ortho)) ^ (b.accept(ortho)));
        }
    }

    public class ByTaxonID implements OrthologueFilter
    {
        int taxonID;

        public ByTaxonID(int taxonID)
        {
            this.taxonID = taxonID;
        }

        public boolean accept(Orthologue ortho)
        {
            try {
                return (taxonID == ortho.getTaxonID());
            }
            catch (Exception e) {
                return false;
            }
        }
    }

    public class ByTaxon implements OrthologueFilter
    {
        Taxon taxon;

        public ByTaxon(Taxon taxon)
        {
            this.taxon = taxon;
        }

        public boolean accept(Orthologue ortho)
        {
            try {
                return (taxon == ortho.getTaxon());
            }
            catch (Exception e) {
                return false;
            }
        }
    }

    public class ByLocusID implements OrthologueFilter
    {
        String locusID;

        public ByLocusID(String locusID)
        {
            this.locusID = locusID;
        }

        public boolean accept(Orthologue ortho)
        {
            try {
                return (locusID.equals(ortho.getLocusID()));
            }
            catch (Exception e) {
                return false;
            }
        }
    }

    public class ByHomologeneID implements OrthologueFilter
    {
        String homologeneID;

        public ByHomologeneID(String homologeneID)
        {
            this.homologeneID = homologeneID;
        }

        public boolean accept(Orthologue ortho)
        {
            try {
                return (homologeneID.equals(ortho.getHomologeneID()));
            }
            catch (Exception e) {
                return false;
            }
        }
    }

    public class ByAccession implements OrthologueFilter
    {
        String regex;

        public ByAccession(String regex)
        {
            this.regex = regex;
        }

        public boolean accept(Orthologue ortho)
        {
            try {
                return ortho.getAccession().matches(regex);
            }
            catch (Exception e) {
                return false;
            }
        }
    }

    public class ByTitle implements OrthologueFilter
    {
        String regex;

        public ByTitle(String regex)
        {
            this.regex = regex;
        }

        public boolean accept(Orthologue ortho)
        {
            try {
                return ortho.getTitle().matches(regex);
            }
            catch (Exception e) {
                return false;
            }
        }
    }
}

