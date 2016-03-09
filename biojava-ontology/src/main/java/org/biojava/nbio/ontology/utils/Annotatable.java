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


package org.biojava.nbio.ontology.utils;



/**
 * <p>Indicates that an object has an associated annotation.</p>
 *
 * <p>Many BioJava objects will have associated unstructured
 * data. This should be stored in an Annotation instance. However, the
 * BioJava object itself will probably not want to extend the
 * Annotation interface directly, but rather delegate off that
 * functionality to an Annotation property. The Annotatable interface
 * indicates that there is an Annotation property. When implementing
 * Annotatable, you should always create a protected or private field
 * containing an instance of ChangeForwarder, and register it as a
 * ChangeListener with the associated Annotation delegate
 * instance.</p>
 *
 * <pre>
 * public class Foo extends AbstractChangeable implements Annotatable {
 *   private Annotation ann; // the associated annotation delegate
 *   protected ChangeForwarder annFor; // the event forwarder
 *
 *   public Foo() {
 *     // make the ann delegate
 *     ann = new SimpleAnnotation();
 *     // construct the forwarder so that it emits Annotatable.ANNOTATION ChangeEvents
 *     // for the Annotation.PROPERTY events it will listen for
 *     annFor = new ChangeForwarder.Retyper(this, getChangesupport( Annotatable.ANNOTATION ),
 *                                          Annotatable.ANNOTATION );
 *     // connect the forwarder so it listens for Annotation.PROPERTY events
 *     ann.addChangeListener( annFor, Annotation.PROPERTY );
 *   }
 *
 *   public Annotation getAnnotation() {
 *     return ann;
 *   }
 * }
 * </pre>
 * Check if BioJava classes and interfaces extend Annotatable. This
 * will tell  you if you should look for associated annotation.
 *
 *  If an object implements Annotatable, it may well propagate
 * ChangeEvent notifications from the associated Annotation. You may
 * need to track these to maintain the state of your applications.
 *
 * Be careful to hook up the appropriate event forwarders.
 *
 * The getAnnotation() method can be implemented lazily
 * (instantiate the Annotation instance and event forwarders when the first
 * request comes in). It can also be implemented by returning throw-away
 * immutable Annotation instances that are built from scratch each time.
 * @author  Matthew Pocock
 * @author <a href="mailto:kdj@sanger.ac.uk">Keith James</a> (docs).
 * @author  Kalle Näslund (docs)
 * @see org.biojavax.RichAnnotatable
 * @since 1.0
 */
public interface Annotatable  {


	/**
	 * Should return the associated annotation object.
	 *
	 * @return an Annotation object, never null
	 */
	Annotation getAnnotation();


}
