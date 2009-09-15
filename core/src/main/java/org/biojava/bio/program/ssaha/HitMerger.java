package org.biojava.bio.program.ssaha;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * A listener that merges overlapping hits and culls all hits under a given
 * length.
 *
 * @author Matthew Pocock
 */
public class HitMerger implements SearchListener {
  private List hitList;
  private int minLength;
  private SearchListener delegate;
  
  /**
   * Build a new HitMerger that will pass events on to a delegate.
   *
   * @param delegate  the SearchListener to inform of all merged and
   *                  filtered hits
   * @param minLength the minimum length a hit must reach to be passed on
   */
  public HitMerger(SearchListener delegate, int minLength) {
    this.hitList = new ArrayList();
    this.delegate = delegate;
    this.minLength = minLength;
  }
  
  public void startSearch(String seqID) {
    hitList.clear();
    delegate.startSearch(seqID);
  }
  
  public void hit(
    int hitID,
    int queryOffset,
    int hitOffset,
    int hitLength
  ) {
    Hit hit = new Hit(hitID, queryOffset, hitOffset, hitLength);
    hitList.add(hit);
  }

  public void endSearch(String seqID) {
    Collections.sort(hitList);

   OUTER:
    for(int i = 0 ; i < hitList.size(); i++) {
      Hit low = (Hit) hitList.get(i);
      int ubound = low.queryOffset + low.hitOffset + low.hitLength * 2;

      for(int j = i+1; j < hitList.size(); j++) {
        Hit high = (Hit) hitList.get(j);

        if(ubound < high.queryOffset + high.hitOffset) {
          break;
        }

        if(doOverlap(low, high)) {
          hitList.set(j, merge(low, high));
          continue OUTER;
        }
      }

      if(low.hitLength >= minLength) {
        delegate.hit(low.hitID, low.queryOffset, low.hitOffset, low.hitLength);
      }
    }

    delegate.endSearch(seqID);
  }

  private Hit merge(Hit a, Hit b) {
    int qo = Math.min(a.queryOffset, b.queryOffset);
    int ho = Math.min(a.hitOffset, b.hitOffset);
    int len = Math.max(
      a.queryOffset + a.hitLength,
      b.queryOffset + b.hitLength ) -
      qo;

    return new Hit(a.hitID, qo, ho, len);
  }

  private boolean doOverlap(Hit a, Hit b) {
    return
      // same sequence
      (a.hitID == b.hitID) &&
      // overlap in query offset space
      ( !(a.queryOffset + a.hitLength < b.queryOffset ||
          a.queryOffset > b.queryOffset + b.hitLength ) ) && 
      // on the same diagonal
      ( a.queryOffset - b.queryOffset == a.hitOffset - b.hitOffset);
  }
  
  private static class Hit
  implements Comparable {
    public int hitID;
    public int queryOffset;
    public int hitOffset;
    public int hitLength;
    
    public Hit(
      int hitID,
      int queryOffset,
      int hitOffset,
      int hitLength
    ) {
      this.hitID = hitID;
      this.queryOffset = queryOffset;
      this.hitOffset = hitOffset;
      this.hitLength = hitLength;
    }
    
    public boolean equals(Object o) {
      if(o instanceof Hit) {
        Hit r = (Hit) o;
        return
          hitID == r.hitID &&
          queryOffset == r.queryOffset &&
          hitOffset == r.hitOffset &&
          hitLength == r.hitLength;
      }
      return false;
    }
    
    public int compareTo(Object o) {
      Hit r = (Hit) o;
      
      if(hitID > r.hitID) {
        return 1;
      } else if(hitID < r.hitID) {
        return -1;
      }
      
      int relDist = queryOffset + hitOffset - (r.queryOffset + r.hitOffset);
      
      if(relDist > 0) {
        return 1;
      } else if(relDist < 0) {
        return -1;
      } else if(hitOffset > r.hitOffset) {
        return 1;
      } else if(hitOffset < r.hitOffset) {
        return -1;
      } else if(hitLength > r.hitLength) {
        return 1;
      } else if(hitLength < r.hitLength) {
        return -1;
      } else {
        return 0;
      }
    }
    
    public String toString() {
      return hitID + " " + queryOffset + " " + hitOffset + " " + hitLength;
    }
  }
}
