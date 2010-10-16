/*
 * BioJava development code
 * 
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 * 
 * http://www.gnu.org/copyleft/lesser.html
 * 
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 * 
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 * 
 * http://www.biojava.org
 */

package org.biojava.utils;

import java.io.PrintStream;
import java.util.EventListener;

/**
 *  Interface for objects which listen to ChangeEvents. 
 *
 * @author     Thomas Down
 * @author     Matthew Pocock
 * @author     Keith James
 * @since      1.1 
 */

public interface ChangeListener extends EventListener {
  /**
   * Convenience implementation which vetoes every change of which it is
   * notified. You could add this to an object directly to stop it changing
   * in any way, or alternatively you could add it for a specific ChangeType
   * to stop that single item from altering.
   */
  final static ChangeListener ALWAYS_VETO = new AlwaysVetoListener();

  /**
   * Convenience implementation that echoes all events to out.
   */
  final static ChangeListener LOG_TO_OUT = new LoggingListener(System.out);
   
  /**
   * <p>
   * Called before a change takes place.
   * </p>
   *
   * <p>
   * This is your chance to stop the change by throwing a ChangeVetoException.
   * This method does not indicate that the change will definitely take place,
   * so it is not recomended that you take any positive action within this
   * handler.
   * </p>
   *
   * @param  cev                      
   * An event encapsulating the change which is about 
   * to take place.
   * @exception  ChangeVetoException  Description of Exception 
   * @throws                          
   * ChangeVetoException if the receiver does not wish 
   * this change to occur at this
   * time.
   */

  void preChange(ChangeEvent cev) throws ChangeVetoException;


  /**
   * <p>
   * Called when a change has just taken place.
   * </p>
   *
   * <p>
   * This method is the place to perform any behavior in response to the
   * change event.
   * </p>
   *
   * @param  cev  
   * An event encapsulating the change which has 
   * occured.
   */

  void postChange(ChangeEvent cev);

  /**
   *  An implementation that always vetoes everything. 
   *
   * @author     Thomas Down
   * @author     Matthew Pocock
   * @since      1.1 
   */

  static class AlwaysVetoListener implements ChangeListener {

    /**
     *  Private constructor.
     */
    protected AlwaysVetoListener() {
    }

    public void preChange(ChangeEvent cev) throws ChangeVetoException {
      throw new ChangeVetoException(
        cev,
        "This object has been locked"
      );
    }

    /**
     * @throws AssertionFailure if this is called, as preChange should have
     *   vetoed any change already
     */
    public void postChange(ChangeEvent cev) {
      throw new AssertionFailure(
        new ChangeVetoException(
          cev,
          "This object has been locked" ));
    }
  }
  
  /**
   * A listener that remembers the ChangeEvent of the last change. Mostly for
   * debugging.
   * @author Mark Schreiber
   * @since 1.5
   */
   
  public class ChangeEventRecorder extends ChangeAdapter{
      private ChangeEvent event;
      
      public ChangeEvent getEvent(){return this.event;}
      
      public void preChange(ChangeEvent cev) throws ChangeVetoException {
          event = cev;
      }
  }
  
  /**
   * A listener that writes information about the event stream to a PrintStream.
   * 
   * @author Matthew Pocock
   * @since 1.1
   */
  public class LoggingListener implements ChangeListener {
    private PrintStream out;
    private String prefix;
    
    /**
     * Create a LoggingListener that will log all events to 'out'.
     *
     * @param out the PrintStream to log events to
     */
    public LoggingListener(PrintStream out) {
      this.out = out;
      this.prefix = null;
    }
    
    /**
     * Create a LoggingListener that will log all events to 'out' with a prefix.
     *
     * @param out the PrintStream to log events to
     * @param prefix the prefix to attach to each line of the log
     */
    public LoggingListener(PrintStream out, String prefix) {
      this.out = out;
      this.prefix = prefix;
    }
    
    public void preChange(ChangeEvent cev) throws ChangeVetoException {
      if(prefix != null) {
        out.print(prefix);
      }
      out.println("preChange for event " + cev);
    }
    
    public void postChange(ChangeEvent cev) {
      if(prefix != null) {
        out.print(prefix);
      }
      out.println("postChange for event " + cev);
    }
  }
}

