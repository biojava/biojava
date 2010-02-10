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

import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.biojava.bio.seq.ByLocationMinMaxFeatureComparator;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;


/**
 * A SequenceRenderer that renders a set of Features that match a FeatureFilter in such a way that
 * they do not overlap in the display.
 *
 * @author Mark Southern
 * @since 1.5
 */
public abstract class AbstractPeptideDigestRenderer extends MultiLineRenderer
{
    public static final ChangeType DIGEST = new ChangeType("The peptide digest has changed",
        "org.biojava.bio.gui.sequence.AbstractPeptideDigestRenderer", "DIGEST",
        SequenceRenderContext.REPAINT
    );
    public static final String LANE = "Lane";
    private FeatureSource source;
    private FeatureFilter digestFilter;
    private Map laneMap = new HashMap();
    private int laneCount = 0;
    private int distanceBetween = 0;
    
    public AbstractPeptideDigestRenderer()
    {
    	super();
    }

    public AbstractPeptideDigestRenderer(FeatureSource source)
    {
    	this();
        setFeatureSource(source);
    }

    public AbstractPeptideDigestRenderer(FeatureSource source, FeatureFilter filter)
    {
        this(source);
        setFilter(filter);
    }
    
    public AbstractPeptideDigestRenderer(FeatureSource source, FeatureFilter filter, int distanceBetweenFeatures)
    {
        this(source,filter);
        setDistanceBetweenFeatures(distanceBetweenFeatures);
    }
    
    public void setFeatureSource(FeatureSource source)
    {
    	this.source = source;
    }
    
    public FeatureSource getFeatureSource()
    {
    	return source;
    }

    public FeatureFilter getFilter()
    {
        return digestFilter;
    }

    public void setFilter(FeatureFilter filter)
    {
        digestFilter = filter;
    }

    /*
     * Sets the space between rendered features. Increase for greater visibility. 
     */
    public void setDistanceBetweenFeatures(int d)
    {
        distanceBetween = d;
    }

    public int getDistanceBetweenFeatures()
    {
        return distanceBetween;
    }

    public void sortPeptidesIntoLanes() throws ChangeVetoException
    {
        if (hasListeners(DIGEST))
        {
            ChangeSupport cs = getChangeSupport(SequenceRenderContext.REPAINT);

            synchronized (cs)
            {
                ChangeEvent ce = new ChangeEvent(this, DIGEST);
                cs.firePreChangeEvent(ce);
                doSortPeptides();
                doRefreshRenderers();
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            doSortPeptides();
            doRefreshRenderers();
        }
    }

    protected void doRefreshRenderers() throws ChangeVetoException
    {
        super.clearRenderers();
        // resort peptide features into new lanes
        for (int j = 1; j <= laneCount; j++)
        {
            //logger.debug("Adding renderers for lane " + j);
            FeatureFilter ffilt = new FeatureFilter.And(getFilter(), new LaneFeatureFilter(j));
            FeatureBlockSequenceRenderer block = new FeatureBlockSequenceRenderer();
            block.setFeatureRenderer(createRenderer(j));
            PaddingRenderer pad = new PaddingRenderer();
            pad.setPadding(1);
            pad.setRenderer(new FilteringRenderer(block, ffilt, true));
            addRenderer(pad);
        }
    }
    
    /*
     * Method used to return the given FeatureRenderer for a given lane in the display.
     */
    public abstract FeatureRenderer createRenderer(int lane);
    
    protected void doSortPeptides()
    {
        // clear existing stored features
        laneMap.clear();

        //logger.debug("Feature Filter = " + getFilter());
        FeatureHolder fh = source.getFeatureHolder().filter(getFilter());
        List ranges = new LinkedList();

        for (Iterator i = fh.features(); i.hasNext();)
        {
            ranges.add(( Feature ) i.next());
        }

        Collections.sort(ranges, new ByLocationMinMaxFeatureComparator());

        Integer lane_id = new Integer(1);
        int i = 0;
        int pos = 0;

        while (ranges.size() > 0)
        {
            /*//logger.info("i=" + i + "\tpos=" + pos + "\tlane_id=" + lane_id + "\tsize=" +
               ranges.size()
               );
             */
            Feature f = ( Feature ) ranges.get(i);

            //logger.info("\tloc=" + f.getLocation());
            if (f.getLocation().getMin() > pos)
            {
                //logger.info("Adding location " + f.getLocation() + " to lane " + lane_id);
                pos = f.getLocation().getMax() + distanceBetween; // +1 so there are no adjoining peptides on a track

                // we can make distanceBetween 0 if we can differenciate between adjoining Features with the specific SequenceRenderer implementation
                ranges.remove(i);
                laneMap.put(f, lane_id);
            } else
            {
                i++;
            }

            if (i >= ranges.size())
            {
                i = 0;
                pos = 0;
                lane_id = new Integer(lane_id.intValue() + 1);
                //logger.info("RESET:\t" + i + "\t" + pos + "\t" + lane_id + "\t" + ranges.size());
            }
        }
        laneCount = lane_id.intValue();
    }

    private class LaneFeatureFilter implements FeatureFilter
    {
        private int lane;

        public LaneFeatureFilter(int lane)
        {
            this.lane = lane;
        }

        public boolean accept(Feature f)
        {
            Integer i = ( Integer ) laneMap.get(f);
            return (i != null) && (i.intValue() == lane);
        }
    }
}
