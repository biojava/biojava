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

package org.biojava.bio.program.tagvalue;

import java.util.HashMap;
import java.util.Map;

import org.biojava.utils.ParserException;

/**
 * This class implements a state machine for parsing events from
 * the Parser class.
 * </p>
 * <p>
 * Each State can be specified to deliver events
 * to a particular TagValueListener.
 * </p>
 * <p>
 * Transitions can be specified to occur between States when
 * specific events are encountered. These events can be
 * Tags delivered by startTag as well as the endTag/endRecord
 * events.  Events that result in exit from the current State
 * can be specified to be notifiable to the State listener.
 * </p>
 * <p>
 * In addition, tables of transitions can be defined and specified
 * as the fallback when the a corresponding Transition cannot be
 * found for the tag.  This is useful for specifying the destination
 * States that are globally applicable for particular tags.  As the
 * fallbacks can be chained, you can end up with a hierarchy of
 * Transition Tables, some being State-specific, others applicable to
 * groups of States and finally one being global.
 * </p>
 * @author David Huen
 */
public class StateMachine
	implements TagValueWrapper
{
	protected static final String START_RECORD_TAG = "__START_RECORD_TAG__";
	protected static final String END_RECORD_TAG = "__END_RECORD__";
	protected static final String END_TAG = "__END_TAG__";
	protected static final String MAGICAL_STATE = "__MAGICAL__";

	protected TagValueListener delegate = null;

	/**
	 * Interface implemented by State listeners that
	 * want notification when a transition leaves the State.
	 */
	public interface ExitNotification
	{
		public void notifyExit() throws ParserException;
	}

	/**
	 * Interface for a State within this StateMachine
	 */
	public interface State
	{
    	public String getLabel();
	    public TagValueListener getListener();
//    	public void setListener(TagValueListener listener);
//	    public void setTransition(Object tag, State destination);
//    	public void setDefaultTransition(State defaultDestination);
	    public void transit(Object tag) throws ParserException;
	}

	/**
	 * class to represent a State Transition
	 */
	public class Transition
	{
		/**
		 * the terminus of this Transition
		 */
		public State destination = null;

		/**
		 * should the State listener be notified
		 * that this transition is leaving the state.
		 */
		public boolean notifyOnExit = false;

		private Transition(State destination, boolean notifyOnExit)
		{
			this.destination = destination;
			this.notifyOnExit = notifyOnExit;
		}
	}

	/**
	 * Table of Transition destination States
	 * and their corresponding Tags.
	 * <p>
	 * Note that you can chain a series of these
	 * Transition tables and the lookup will
	 * proceed along the chain until it succeeds.
	 */
	public class TransitionTable
	{
		private Map table = new HashMap();
		private TransitionTable fallback = null;

		protected TransitionTable() {}

		/**
		 * set a Transition within this TransitionTable (2-argument form)
		 */
		public void put(Object tag, Transition transition)
			throws ParserException
		{
			if (table.containsKey(tag))
				throw new ParserException("duplicate tag " + tag);

			table.put(tag, transition);
		}

		/**
		 * set a Transition within this TransitionTable (3-argument form)
		 */
		public void setTransition(Object tag, State destination, boolean notifyOnExit)
			throws ParserException
		{
			put(tag, new Transition(destination, notifyOnExit));
		}

		/**
		 * get the Transition associated with the specified tag.
		 */
		public Transition get(Object tag)
		{
			Transition transition =  (Transition) table.get(tag);

			if ((transition == null) && (fallback != null)) {
				// do a fallback lookup
				transition = fallback.get(tag);
			}

			return transition;
		}

		/**
		 * set the specified TransitionTable to be looked
		 * looked up if the Transition cannot be found in
		 * this one.
		 */
		public void setFallback(TransitionTable fallback)
		{
			this.fallback = fallback;
		}
	}


	/**
	 * Implementation of a State in a state machine
	 */
	public class BasicState
		implements State
	{
		private String label;
		private TransitionTable transitions;
		private	TransitionTable defaultTransitions;
		private TagValueListener listener = null;

		/**
		 * This is the default constructor
		 */
		public BasicState(String label)
		{
			this.label = label;
			this.transitions = new TransitionTable();
		}

		/**
		 * when this constructor is used, a fixed listener
		 * is used with this state.  When setListener is called
		 * on this class, it is the listener of this fixed
		 * listener that is changed.
		 */
		public BasicState(String label, TagValueListener listener)
		{
			this.label = label;
			this.listener = listener;
			this.transitions = new TransitionTable();
		}

		/**
		 * return the label of this class.
		 */
		public String getLabel() { return label; }

		/**
		 * return the TagValueListener assigned to this State.
		 */
		public TagValueListener getListener() { return listener; }

		/**
		 * set a TagValueListener for this State.
		 */
		public void setListener(TagValueListener listener)
		{
			this.listener = listener;
		}

		/**
		 * set a Transition for this State
		 */
		public void setTransition(Object tag, State destination, boolean notifyOnExit)
			throws ParserException
		{
			transitions.setTransition(tag, destination, notifyOnExit);
		}

		/**
		 * set a Transition for this State setting notifyOnExit to false.
		 */
		public void setTransition(Object tag, State destination)
			throws ParserException
		{
			setTransition(tag, destination, false);
		}

		/**
		 * retrieve the TransitionTable for this State.
		 */
		public TransitionTable getTransitionTable()
		{
			return transitions;
		}

		/**
		 * specify fallback TransitionTable for this State
		 */
		public void setDefaultTransitions(TransitionTable defaultTransitions)
		{
			this.defaultTransitions = defaultTransitions;
		}

		/**
		 * Find the destination State when the specified tag
		 * is encountered.  If the destination is successfully
		 * obtained, the listener is also notified of exit if
		 * required.
		 */
		public void transit(Object tag)
			throws ParserException
		{

			Transition transition = transitions.get(tag);

			if (transition == null) {
				// explicit transition unavailable.
				// Use fallback default.
				if (defaultTransitions == null) 
					throw new ParserException("no transition available from " + statePointer.getLabel() + " for tag '" + tag + "'");

				transition = defaultTransitions.get(tag);

				if (transition == null)
					throw new ParserException("no transition available from " + statePointer.getLabel() + " for tag '" + tag + "'");
			}

			State nextState = transition.destination;

			if (transition.notifyOnExit) {
//				System.out.println("notifyOnExit being checked. this is " + this.getLabel() + " next is " + nextState.getLabel());
				if (nextState != this) {
					if (listener instanceof ExitNotification) {
						((ExitNotification) listener).notifyExit();
					}
				}
			}

			statePointer = nextState;
		}
	}

	/**
	 * a basic listener for a State.  It forwards all events to the
	 * delegate for the StateMachine.  Extend to implement listeners
	 * for specific states.
	 */
	public class SimpleStateListener
		implements TagValueListener
	{
		private boolean exceptionOnNullDelegate = true;

		/**
		 * determines if an exception is thrown when an event
		 * arrives without the delegate being set.  Default
		 * is that a ParserException is thrown.
		 i*/
		public void setExceptionOnNullDelegate(boolean throwException)
		{ exceptionOnNullDelegate = throwException; }

		public void startTag(Object tag)
			throws ParserException 
		{ 
			if (delegate != null)
				delegate.startTag(tag);
			else if (exceptionOnNullDelegate)
				throw new ParserException("event arrived without a delegate being specified");
		}

		public void endTag() 
			throws ParserException 
		{
			if (delegate != null)
		 		delegate.endTag();
			else if (exceptionOnNullDelegate)
				throw new ParserException("event arrived without a delegate being specified");
		}

		public void startRecord() 
			throws ParserException 
		{ 
			if (delegate != null)
				delegate.startRecord();
			else if (exceptionOnNullDelegate)
				throw new ParserException("event arrived without a delegate being specified"); 
		}

		public void endRecord() 
			throws ParserException 
		{ 
			if (delegate != null)
				delegate.endRecord();
			else if (exceptionOnNullDelegate)
				throw new ParserException("event arrived without a delegate being specified"); 
		}

		public void value(TagValueContext ctxt, java.lang.Object value)
			throws ParserException
		{ 
			if (delegate != null)
				delegate.value(ctxt, value);
			else if (exceptionOnNullDelegate) 
				throw new ParserException("event arrived without a delegate being specified");
		}
	}


	/**
	 * current state
	 */
	private State statePointer;

	/**
	 * mapping between state labels and states
	 */
	private Map states = new HashMap();
	private BasicState magicalState;

	public StateMachine()
	{
		magicalState = createState(MAGICAL_STATE);
		statePointer = magicalState;
	}

	public BasicState getMagicalState() { return magicalState; }

	public BasicState createState(String label)
	{
		if (states.get(label)!= null) {
			throw new IllegalArgumentException("label " + label + " is not unique");
		}

		BasicState newState = new BasicState(label);
		states.put(label, newState);

		return newState;
	}

	public State getState(String label)
	{
		return (State) states.get(label);
	}

	public TransitionTable createTransitionTable()
	{
		return new TransitionTable();
	}

	/**
	 * TagValueWrapper interface
	 */
	public void setDelegate(TagValueListener delegate) { this.delegate = delegate; }
	public TagValueListener getDelegate() { return delegate; }

	/**
	 * TagValueListener interface
	 */
	public void startTag(Object tag)
		throws ParserException
	{
		statePointer.transit(tag);

		TagValueListener listener = statePointer.getListener();
		if (listener != null) listener.startTag(tag);
	}

	public void endTag()
		throws ParserException
	{
		TagValueListener lstnr;
		if ((lstnr = ((BasicState) statePointer).getListener()) != null)
			lstnr.endTag();
		statePointer.transit(END_TAG);
	}

	public void startRecord()
		throws ParserException
	{
		TagValueListener lstnr;
		if ((lstnr = ((BasicState) statePointer).getListener()) != null)
			lstnr.startRecord();
		statePointer.transit(START_RECORD_TAG);
	}

	public void endRecord()
		throws ParserException
	{
		TagValueListener lstnr;
		if ((lstnr = ((BasicState) statePointer).getListener()) != null)
			lstnr.endRecord();
		statePointer.transit(END_RECORD_TAG);
	}

	public void value(TagValueContext ctxt, Object value)
		throws ParserException
	{
		TagValueListener lstnr;
		if ((lstnr = ((BasicState) statePointer).getListener()) != null)
			lstnr.value(ctxt, value);
	}
}
