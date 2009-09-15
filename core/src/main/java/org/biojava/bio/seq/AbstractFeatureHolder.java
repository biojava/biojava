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

package org.biojava.bio.seq;

import java.util.Iterator;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeVetoException;

/**
 * An abstract implementation of FeatureHolder.
 *
 * This provides the filter method, but who wants to code that more than
 * once? It also has support for the ChangeEvents.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */
public abstract class AbstractFeatureHolder
  extends
    AbstractChangeable
  implements
    FeatureHolder
{
    public FeatureHolder filter(FeatureFilter filter) {
        boolean recurse = !FilterUtils.areProperSubset(filter, FeatureFilter.top_level);
        return filter(filter, recurse);
    }

  public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
    SimpleFeatureHolder res = new SimpleFeatureHolder();
    for(Iterator f = features(); f.hasNext();) {
      Feature feat = (Feature) f.next();
      if(ff.accept(feat)) {
        try {
          res.addFeature(feat);
        } catch (ChangeVetoException cve) {
          throw new BioError(
            "Assertion failed: Couldn't add a feature to my new FeatureHolder"
          );
        }
      }
      if(recurse) {
        FeatureHolder r = feat.filter(ff, recurse);
        for(Iterator rf = r.features(); rf.hasNext();) {
          try {
            res.addFeature((Feature) rf.next());
          } catch (ChangeVetoException cve) {
            throw new BioError(
              "Assertion failure: Should be able to manipulate this FeatureHolder", cve
            );
          }
        }
      }
    }
    return res;
  }

  public Feature createFeature(Feature.Template temp)
  throws BioException, ChangeVetoException {
    throw new ChangeVetoException(
      "This FeatureHolder does not support creation of new Features."
    );
  }

  public void removeFeature(Feature f)
  throws ChangeVetoException, BioException {
    throw new ChangeVetoException(
      "This FeatureHolder does not support removal of Features."
    );
  }
}
