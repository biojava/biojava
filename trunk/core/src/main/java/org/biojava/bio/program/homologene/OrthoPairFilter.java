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


public interface OrthoPairFilter
{
    public boolean accept(OrthoPair pair);


    public final static class AcceptAll implements OrthoPairFilter
    {
        public boolean accept(OrthoPair pair) { return true; }
    }

    public final static class Not implements OrthoPairFilter
    {
        OrthoPairFilter a;

        public Not(OrthoPairFilter a)
        {
            this.a = a;
        }

        public boolean accept(OrthoPair pair)
        {
            return !a.accept(pair);
        }
    }

    public final static class Or implements OrthoPairFilter
    {
        OrthoPairFilter a;
        OrthoPairFilter b;

        public Or(OrthoPairFilter a, OrthoPairFilter b)
        {
            this.a = a;
            this.b = b;
        }

        public boolean accept(OrthoPair pair)
        {
            return ((a.accept(pair)) || (b.accept(pair)));
        }
    }

    public final static class And implements OrthoPairFilter
    {
        OrthoPairFilter a;
        OrthoPairFilter b;

        public And(OrthoPairFilter a, OrthoPairFilter b)
        {
            this.a = a;
            this.b = b;
        }

        public boolean accept(OrthoPair pair)
        {
            return ((a.accept(pair)) && (b.accept(pair)));
        }
    }

    public final static class Xor implements OrthoPairFilter
    {
        OrthoPairFilter a;
        OrthoPairFilter b;

        public Xor(OrthoPairFilter a, OrthoPairFilter b)
        {
            this.a = a;
            this.b = b;
        }

        public boolean accept(OrthoPair pair)
        {
            return ((a.accept(pair)) ^ (b.accept(pair)));
        }
    }

    public class ByMinIdentity implements OrthoPairFilter
    {
        double lowLimit;

        public ByMinIdentity(double lowLimit)
        {
            this.lowLimit = lowLimit;
        }

        public boolean accept(OrthoPair pair)
        {
            // it is possible that an OrthoPair
            // may not have the percentIdentity
            // defined, e.g. curated orthologies
            // the  method will return false then.

            // check that the SimilarityType is acceptable
            SimilarityType pairType = pair.getSimilarity();

            if ( (pairType == SimilarityType.TWIN)
                   || (pairType == SimilarityType.MULTIPLE) ) {

                return (pair.getPercentIdentity() >= lowLimit);
            }
            else
                return false;
        }
    }

    public class ByMaxIdentity implements OrthoPairFilter
    {
        double hiLimit;

        public ByMaxIdentity(double hiLimit)
        {
            this.hiLimit = hiLimit;
        }

        public boolean accept(OrthoPair pair)
        {
            // it is possible that an OrthoPair
            // may not have the percentIdentity
            // defined, e.g. curated orthologies
            // the  method will return false then.

            // check that the SimilarityType is acceptable
            SimilarityType pairType = pair.getSimilarity();

            if ( (pairType == SimilarityType.TWIN)
                   || (pairType == SimilarityType.MULTIPLE) ) {

                return (pair.getPercentIdentity() < hiLimit);
            }
            else
                return false;
        }
    }

    public class BySimilarityType implements OrthoPairFilter
    {
        SimilarityType type;

        public BySimilarityType(SimilarityType type)
        {
            this.type = type;
        }

        public boolean accept(OrthoPair pair)
        {
            return (pair.getSimilarity() == type);
        }
    }

    public class ByRef implements OrthoPairFilter
    {
        String regex;

        public ByRef(String regex)
        {
            this.regex = regex;
        }

        public boolean accept(OrthoPair pair)
        {
            if (pair.getSimilarity() == SimilarityType.CURATED) {
                return pair.getRef().matches(regex);
            }
            else 
                return false;
        }
    }

}
