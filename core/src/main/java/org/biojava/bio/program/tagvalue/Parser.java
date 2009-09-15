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

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.utils.ParserException;

/**
 * <p>
 * Encapsulate the parsing of lines from a buffered reader into tag-value
 * events.
 * </p>
 *
 * <p>
 * Scripts will usually construct a parser object, a BufferedReader, a
 * TagValueParser and TagValueListener, and then set up a loop that processes
 * each record in the reader by calling Parser.read() until it returns false.
 * </p>
 *
 * @since 1.2
 * @author Matthew Pocock
 */
public class Parser {
  public boolean read(
    BufferedReader reader,
    TagValueParser parser,
    TagValueListener listener
  ) throws
    IOException,
    ParserException
  {
    final List stack = new ArrayList();
    push(stack, new Frame(parser, listener, null));
    
    Context ctxt = new Context();
    
    listener.startRecord();
    
    for(
      Object line = reader.readLine();
      line != null;
      line = reader.readLine()
    ) {
      // Find the deepest stack-frame with the same key.
      // If this is not the deepest one, unwind the stack to that point
      Frame frame = null;
      TagValue tv = null;
      for(Iterator fi = stack.iterator(); fi.hasNext(); ) {
        frame = (Frame) fi.next();
        
        tv = frame.parser.parse(line);
        
        // end of record. Is it the last in the file?
        if(tv == null) {
          // scan for eof after whitespace
          boolean eof = false;
          while(true) {
            reader.mark(1);
            int c = reader.read();
            if(c == -1) {
              eof = true;
              break;
            }
            
            if(Character.isWhitespace((char) c)) {
              continue;
            }
            
            reader.reset();
            break;
          }
          
          // now unwind stack
          do {
            Frame top = (Frame) pop(stack);
            top.listener.endTag();
            top.listener.endRecord();
          } while(!stack.isEmpty());
          return !eof;
        }
        
        // not a continuation of a previous tag - unwrap stack
        if(
          tv.isNewTag() ||
          (tv.getTag() != null && !tv.getTag().equals(frame.tag))
        ) {
          // remove all stack frames which have been obsoleted by this tag
          Frame top;
          for(top = (Frame) pop(stack); top != frame; top = (Frame) pop(stack)) {
            top.listener.endTag();
            top.listener.endRecord();
          }
          if(top.tag != null) {
            top.listener.endTag();
          }
          
          // handle current stack frame by starting a tag
          push(stack, new Frame(top.parser, top.listener, tv.getTag()));
          top.listener.startTag(tv.getTag());
          break;
        }
        
        line = tv.getValue();
      }
      
      // process a value and handle potentially pushing a new stack frame
      while(true) {
        // pass in value and see if it requests a new stack frame
        ctxt.flush();
        frame.listener.value(ctxt, tv.getValue());
        if(!ctxt.isDirty()) {
          break;
        }

        // push a new stack frame
        tv = ctxt.parser.parse(tv.getValue());
        frame = new Frame(ctxt.parser, ctxt.listener, tv.getTag());
        push(stack, frame);
        ctxt.listener.startRecord();
        ctxt.listener.startTag(tv.getTag());
        
        // we must loop arround in case the new frame wants to immediately push a
        // new stack frame
      }
    }

    throw new IOException("Premature end of stream or missing end tag");
  }
  
  private void push(List stack, Object o) {
    stack.add(o);
  }
  
  private Object pop(List stack) {
    return stack.remove(stack.size() - 1);
  }
  
  private static class Frame {
    public final TagValueParser parser;
    public final TagValueListener listener;
    public final Object tag;
    
    public Frame(
      TagValueParser parser,
      TagValueListener listener,
      Object tag
    ) {
      this.parser = parser;
      this.listener = listener;
      this.tag = tag;
    }
  }
  
  private static class Context
    implements
      TagValueContext
  {
    public TagValueParser parser;
    public TagValueListener listener;
    
    public void pushParser(
      TagValueParser subParser,
      TagValueListener listener
    ) {
      this.parser = subParser;
      this.listener = listener;
    }
    
    public void flush() {
      this.parser = null;
      this.listener = null;
    }
    
    public boolean isDirty() {
      return parser != null;
    }
  }
}
