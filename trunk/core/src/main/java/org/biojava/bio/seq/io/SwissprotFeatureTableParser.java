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


package org.biojava.bio.seq.io;

import java.util.StringTokenizer;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.symbol.FuzzyLocation;
import org.biojava.bio.symbol.FuzzyPointLocation;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;

/**
 * Simple parser for swissprot feature tables.
 *
 * @author Greg Cox
 * @author Thomas Down
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 */

/*
 * Thomas Down refactored to remove dependancy on FeatureTableParser, since we
 *             weren't actually reusing very much code, and FeatureTableParser
 *             changed quite a bit when it went fully-newio.
 */

class SwissprotFeatureTableParser
{
    private SeqIOListener listener;
    private String featureSource;

    private boolean inFeature = false;
    private Feature.Template featureTemplate;
    private StringBuffer descBuf;

    {
        descBuf = new StringBuffer();
    }

    SwissprotFeatureTableParser(SeqIOListener listener, String source)
    {
        this.listener = listener;
        this.featureSource = source;
    }

    public void startFeature(String type)
        throws BioException
    {
        featureTemplate = new Feature.Template();
        featureTemplate.source = featureSource;
        featureTemplate.type = type;
        descBuf.setLength(0);
        inFeature = true;
    }

        public void featureData(String line)
            throws BioException
        {
                boolean newFeature = false;
                // Check if there is a location section.
                if(line.charAt(5) != ' ')
                {
                        StringTokenizer tokens = new StringTokenizer(line);
                        featureTemplate.location = getLocation(tokens);

                        if(line.length() >= 20)
                        {
                                line = line.substring(20);
                        }
                        else
                        {
                                line = "";
                        }
                        newFeature = true;
                }

                if(newFeature == true)
                {
                        descBuf.setLength(0);
                }
                descBuf.append(" " + line.trim());
                newFeature = false;
        }

        public void endFeature()
            throws BioException
        {
            if (descBuf.length() > 0) {
                featureTemplate.annotation = new SimpleAnnotation();
                try {
                    featureTemplate.annotation.setProperty(SwissprotProcessor.PROPERTY_SWISSPROT_FEATUREATTRIBUTE, descBuf.substring(0));
                } catch (ChangeVetoException ex) {
                    throw new BioException("Couldn't alter annotation",ex);
                }
            } else {
                featureTemplate.annotation = Annotation.EMPTY_ANNOTATION;
            }

            listener.startFeature(featureTemplate);
            listener.endFeature();

            inFeature = false;
        }

        public boolean inFeature()
        {
            return inFeature;
        }

        /**
         * Returns the next location contained in theTokens
         *
         * @exception bioException Thrown if a non-location is first in theTokens
         * @param theTokens The tokens to process
         * @return The location at the front of theTokens
         */
        private Location getLocation(StringTokenizer theTokens)
                throws BioException
        {
                Index startIndex = this.getIndex(theTokens);
                Index endIndex = this.getIndex(theTokens);

                Location theLocation;
                if(startIndex.isFuzzy() || endIndex.isFuzzy())
                {
                        // This handles all locations with one of the following points:
                        // 		<nn
                        //		>nn
                        //		?nn
                        //		?
                        // All these get resolved into a fuzzy location, though with
                        // different combinations of start and end points.  The getIndex
                        // method does some preliminary work to find mins and maxes
                        if(startIndex.equals(endIndex))
                        {
                                theLocation = new FuzzyPointLocation(startIndex.getMinValue(),
                                                startIndex.getMaxValue(),
                                                FuzzyPointLocation.RESOLVE_AVERAGE);
                        }
                        else
                        {
                                theLocation = new FuzzyLocation(startIndex.getMinValue(),
                                        endIndex.getMaxValue(), startIndex.getMaxValue(),
                                        endIndex.getMinValue(), startIndex.isFuzzy(),
                                        endIndex.isFuzzy(), FuzzyLocation.RESOLVE_INNER);
                        }
                }
                else if(startIndex.equals(endIndex))
                {
                        // Fuzzy point locations were peeled off as fuzzy locations.  The
                        // point locations handled here are concrete point locations
                    theLocation = new PointLocation(startIndex.getMinValue());
                }
            else
                {
                    theLocation = new RangeLocation(startIndex.getMinValue(), endIndex.getMinValue());
                }

                return theLocation;
        }

        /**
         * Returns the Integer value of the next token and its fuzzyness
         *
         * @exception BioException Thrown if a non-number token is passed in
         * (fuzzy locations are handled)
         * @param theTokens The tokens to be processed
         * @return Index
         */
        private Index getIndex(StringTokenizer theTokens)
                throws BioException
        {
                String returnIndex = theTokens.nextToken();
                int minValue;
                int maxValue;
                boolean isFuzzy;
                if(returnIndex.indexOf('<') != -1)
                {
                        // Peel off cases "<nn"  Returned as MIN_VALUE to nn
                        returnIndex = returnIndex.substring(1);
                        minValue = Integer.MIN_VALUE;
                        maxValue = Integer.parseInt(returnIndex);
                        isFuzzy = true;
                }
                else if(returnIndex.indexOf('>') != -1)
                {
                        // Peel off cases ">nn"  Returned as nn to MAX_VALUE
                        returnIndex = returnIndex.substring(1);
                        maxValue = Integer.MAX_VALUE;
                        minValue = Integer.parseInt(returnIndex);
                        isFuzzy = true;
                }
                else if(returnIndex.indexOf('?') != -1)
                {
                        // Two cases are handled here; '?' and "?nn"
                        if(returnIndex.length() == 1)
                        {
                                minValue = Integer.MIN_VALUE;
                                maxValue = Integer.MAX_VALUE;
                                isFuzzy = true;
                        }
                        else
                        {
                                returnIndex = returnIndex.substring(1);
                                minValue = Integer.parseInt(returnIndex);
                                maxValue = minValue;
                                isFuzzy = true;
                        }
                }
                else
                {
                        // Plain vanilla location
                        minValue = Integer.parseInt(returnIndex);
                        maxValue = minValue;
                        isFuzzy = false;
                }

                return(new Index(minValue, maxValue, isFuzzy));
        }

        /**
         * This inner class is a struct to pass back the information contained in
         * a Swissprot point
         */
        private class Index
        {
                protected int mMinValue;
                protected int mMaxValue;
                protected boolean isFuzzy;

                public Index(int theMinValue, int theMaxValue, boolean theFuzzyness)
                {
                        this.mMinValue = theMinValue;
                        this.mMaxValue = theMaxValue;
                        this.isFuzzy = theFuzzyness;
                }

                public int getMaxValue()
                {
                        return this.mMaxValue;
                }

                public int getMinValue()
                {
                        return this.mMinValue;
                }

                public boolean isFuzzy()
                {
                        return this.isFuzzy;
                }

                public boolean equals(Index theIndex)
                {
                        boolean returnValue = false;
                        if((this.mMinValue == theIndex.getMinValue()) &&
                                (this.mMaxValue == theIndex.getMaxValue()) &&
                                (this.isFuzzy == theIndex.isFuzzy()))
                        {
                                returnValue = true;
                        }
                        return returnValue;
                }
        }
}
