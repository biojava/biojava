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

package org.biojava3.core.util;

import java.io.IOException;

/**
 * Simple interface for building XML documents.
 *
 * @author Thomas Down
 * @since 1.3
 */

public interface XMLWriter {
    /**
     * Send raw data to the stream.  Mainly useful for things like DOCTYPE
     * declarations.  Use with care!
     *
     * @param s a string of data to include verbatim in the XML stream
     */
    
    public void printRaw(String s) throws IOException;
    
    /**
     * Open a new namespace-qualified XML tag.
     *
     * @param nsURI A URI for the namespace to use
     * @param localName The name of the tag
     */
    
    public void openTag(String nsURI, String localName) throws IOException;
    
    /**
     * Open a new unqualified XML tag.  This may also be used if you want
     * to do namespace management yourself, independantly of the XMLWriter
     *
     * @param name The name of the tag.
     */
    
    public void openTag(String name) throws IOException;
    
    /**
     * Add an attribute to an element.  This will throw an exception if it's not
     * called immediately after an <code>openTag</code> command.
     *
     * @param nsURI A URI for the namespace to use
     * @param localName The name of the attribute
     * @param value The textual value of the attribute
     */
    
    public void attribute(String nsURI, String localName, String value) throws IOException;
    
    /**
     * Add an un-qualified attribute to an element.  This will throw an exception if it's not
     * called immediately after an <code>openTag</code> command.
     *
     * @param qName The name of the attribute to set
     * @param value The textual value of the attribute
     */
    
    public void attribute(String qName, String value) throws IOException;
    
    /**
     * Prints some textual content in an element.
     */
    
    public void print(String data) throws IOException;
    
    /**
     * Prints some textual content, terminated with a newline character.
     */
    
    public void println(String data) throws IOException;
    
    /**
     * Closes an element
     *
     * @param nsURI A URI for the namespace to use
     * @param qName The name of the tag
     */
    
    public void closeTag(String nsURI, String qName) throws IOException;
    
    /**
     * Closes an un-qualified element.
     * 
     * @param name The tag name
     */
    
    public void closeTag(String name) throws IOException;
    
    /**
     * Hints that a namespace is going to be used in a sub-tree.  Use this method
     * to avoid namespaces that are used only in leaf-nodes of a tree being re-defined
     * every time they are used.  The XMLWriter will generally try to use the suggested
     * prefix for this namespace, but there is no guarentee of this.  In particular, if
     * the namespace is already in use, the current prefix will still be used.  Similarly
     * if the suggested prefix has already been used for another namespace, a new one
     * will be auto-generated.
     *
     * @param nsURI The namespace to declare
     * @param prefixHint A suggested prefix-string for this namespace.
     */
     
    public void declareNamespace(String nsURI, String prefixHint) throws IOException;
    
    /**
     * Close this XMLWriter, and it's underlying stream.
     *
     * @since 1.4
     */
    
    public void close() throws IOException;
}
