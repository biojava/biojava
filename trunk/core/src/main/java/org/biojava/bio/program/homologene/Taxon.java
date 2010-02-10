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

import java.util.HashSet;
import java.util.Set;

public interface Taxon
{
    static Set taxa = new HashSet();

    public static Taxon A_THALIANA = new TaxonStub(3702, "Arabidopsis thaliana");
    public static Taxon B_TAURUS = new TaxonStub(9913, "Bos taurus");
    public static Taxon C_ELEGANS = new TaxonStub(6239, "Caenorhabditis elegans");
    public static Taxon D_RERIO = new TaxonStub(7955, "Danio rerio");
    public static Taxon D_MELANOGASTER = new TaxonStub(7227, "Drosophila melanogaster");
    public static Taxon H_SAPIENS = new TaxonStub(9606, "Homo sapiens");
    public static Taxon H_VULGARE = new TaxonStub(4513, "Hordeum vulgare");
    public static Taxon L_ESCULENTUM = new TaxonStub(4081, "Lycopersicon esculentum");
    public static Taxon M_TRUNCULATA = new TaxonStub(3880, "Medicago truncatula");
    public static Taxon M_MUSCULUS = new TaxonStub(10090, "Mus musculus");
    public static Taxon O_SATIVA = new TaxonStub(4530, "Oryz sativa");
    public static Taxon R_NORVEGICUS = new TaxonStub(10116, "Rattus norvegicus");
    public static Taxon S_SCROFA = new TaxonStub(9823, "Sus scrofa");
    public static Taxon T_AESTIVUM = new TaxonStub(4565, "Triticum aestivum");
    public static Taxon X_LAEVIS = new TaxonStub(8355, "Xenopus laevis");
    public static Taxon Z_MAYS = new TaxonStub(4577, "Zea mays");


    /**
     * returns the name of the Taxon
     */
    public String getDescription();

    /**
     * returns the taxon ID
     */
    public int getTaxonID();

    public static class TaxonStub implements Taxon
    {
        int id;
        private String description;

        TaxonStub(int id, String description)
        {
            this.id = id;
            this.description = description;

            taxa.add(this);
        }

        public String getDescription() { return description; }

        public int getTaxonID() { return id; }

        public int hashCode() { return id; }
    }

}

