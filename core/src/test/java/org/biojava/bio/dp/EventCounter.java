package org.biojava.bio.dp;

import java.util.ArrayList;
import java.util.List;

import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeVetoException;

/**
 * Utility class for event handler tests.
 *
 * @author Matthew Pocock
 */
class EventCounter
implements ChangeListener {
  private final String name;
  private final List preList;
  private final List postList;

  public EventCounter(String name) {
    this.name = name;
    this.preList = new ArrayList();
    this.postList = new ArrayList();
  }

  public void preChange(ChangeEvent cev) throws ChangeVetoException {
    preList.add(cev);
  }

  public void postChange(ChangeEvent cev) {
    postList.add(cev);
  }

  public int getPreCounts() {
    return preList.size();
  }

  public int getPostCounts() {
    return postList.size();
  }

  public void zeroCounts() {
    preList.clear();
    postList.clear();
  }

  public String toString() {
    return "EventCounter(" + name + ")\n\t" + preList + "\n\t" + postList;
  }
}
