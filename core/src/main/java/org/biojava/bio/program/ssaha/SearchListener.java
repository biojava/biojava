package org.biojava.bio.program.ssaha;

import java.io.PrintStream;

/**
 * The interface used to inform interested parties that some sequence has
 * been searched and something found.
 * <p>
 * The callbacks will always be called in the order startSearch, hit,
 * endSearch, during which time there may be multiple hit calls. The seqID
 * of startSearch and endSearch will match. After this, a new startSearch
 * may begin. These events will usually originate from the search method of
 * DataStore.
 *
 * @author Matthew Pocock
 */
public interface SearchListener {
  /**
   * Indicates that a sequence is about to be searched against a DataStore.
   *
   * @param seqID  the id of the sequence to be searched
   */
  public void startSearch(String seqID);

  /**
   * Indicates that a sequence has been searched against a DataStore.
   *
   * @param seqID  the id of the sequence to be searched
   */
  public void endSearch(String seqID);
  
  /**
   * There has been a hit between the query sequence and a database
   * sequence.
   *
   * @param hitID  the number of the sequence hit; resolvable by
   *               String id = DataStore.seqNameForID(hitID)
   * @param queryOffset the offset into the query sequence
   * @param hitOffset the offset into the sequence hit in the database
   * @param hitLength the number of symbols hit
   */
  public void hit(
    int hitID,
    int queryOffset,
    int hitOffset,
    int hitLength
  );

  /**
   * A simple wrapper implementation.
   *
   * <p>Extend this and over-ride any of the interface methods to implement
   * SearchListeners that filter hits before passing them on to an
   * underlying listener.</p>
   * You can modify the search events the delegate sees by over-riding any of
   * the SearchListener methods, modify the arguments
   * and then call the method on super with the new arguments.
   * You can drop hits by just not passing them onto the delegate using
   * super.hits().
   * <em>Note:</em> Be sure to maintain the nesting of start/stop search and
   * hit, or you will confuse the delegate.
   *
   * @author Matthew Pocock
   * @since 1.4
   */
  public static abstract class Wrapper
  implements SearchListener {
    private final SearchListener delegate;

    public Wrapper(SearchListener delegate) {
      this.delegate = delegate;
    }

    public void startSearch(String seqID) {
      delegate.startSearch(seqID);
    }

    public void endSearch(String seqID) {
      delegate.endSearch(seqID);
    }

    public void hit(
      int hitID,
      int queryOffset,
      int hitOffset,
      int hitLength
    ) {
      delegate.hit(hitID, queryOffset, hitOffset, hitLength);
    }
  }

  /**
   * A SearchListener that passes events on to two delegate listeners.
   *
   * <p>This allows you to build trees of listeners. This is usefull, for
   * example, when echoing output from different listeners.</p>
   *
   * @author Matthew Pocock
   * @since 1.4
   */
  public static final class Tee
  implements SearchListener {
    private final SearchListener d1;
    private final SearchListener d2;

    public Tee(SearchListener d1, SearchListener d2) {
      if(d1 == null || d2 == null) {
        throw new IllegalArgumentException(
          "Delegates can not be null: " + d1 + " " + d2 );
      }

      this.d1 = d1;
      this.d2 = d2;
    }

    public void startSearch(String seqID) {
      d1.startSearch(seqID);
      d2.startSearch(seqID);
    }

    public void endSearch(String seqID) {
      d1.endSearch(seqID);
      d2.endSearch(seqID);
    }

    public void hit(
      int hitID,
      int queryOffset,
      int hitOffset,
      int hitLength
    ) {
      d1.hit(hitID, queryOffset, hitOffset, hitLength);
      d2.hit(hitID, queryOffset, hitOffset, hitLength);
    }
  }

  /**
   * A simple listener that filters out all hits that are too short.
   *
   * @author Matthew Pocock
   * @since 1.4
   */
  public static final class FilterByLength
  extends Wrapper {
    private final int minLength;

    public FilterByLength(SearchListener delegate, int minLength) {
      super(delegate);
      this.minLength = minLength;
    }

    public void hit(
      int hitID,
      int queryOffset,
      int hitOffset,
      int hitLength
    ) {
      if(hitLength >= minLength) {
        super.hit(hitID, queryOffset, hitOffset, hitLength);
      }
    }
  }

  /**
   * A SearchListener that prints events out to a PrintStream.
   *
   * <p>Use this for debugging purposes.</p>
   *
   * @author Matthew Pocock
   * @since 1.4
   */
  public static final class Echo
  implements SearchListener {
    private final PrintStream out;

    public Echo(PrintStream out) {
      this.out = out;
    }

    public void startSearch(String seqID) {
      out.println("startSearch: " + seqID);
    }

    public void endSearch(String seqID) {
      out.println("endSearch: " + seqID);
    }

    public void hit(
      int hitID,
      int queryOffset,
      int hitOffset,
      int hitLength
    ) {
      out.println(
        "hit." +
        "\thitID: " + hitID +
        "\tqueryOffset: " + queryOffset +
        "\thitOffset: " + hitOffset +
        "\thitLength: " + hitLength );
    }
  }
}
