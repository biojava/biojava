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
package org.biojava3.core.io;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Puts parser chains together. Parser chains have a {@link ThingReader}
 * or {@link ThingEmitter} at the start which output data about a {@link Thing}.
 * An optional number of {@link ThingConverter}s convert that data into other
 * kinds of {@link Thing}. A final {@link ThingReceiver} which is either a
 * {@link ThingWriter} (writes to files) or a {@link ThingBuilder} (builds new
 * {@link Thing}s) receives the data. Depending on whether the receiver is a
 * writer or a builder, the user can iterate over the data as it is parsed, or
 * have it all written at once to file.
 * @author Richard Holland
 * @since 3.0
 */
public class ThingParser<T extends Serializable> implements Iterator<T> {

    private final Logger logger = Logger.getLogger("org.biojava3.core.io.ThingParser");
    private ThingReader thingReader;
    private List<ThingConverter<? extends Serializable>> thingConverters = new ArrayList<ThingConverter<? extends Serializable>>();
    private ThingReceiver thingReceiver;
    private Exception lastException;

    /**
     * Construct a parser that reads from one source and writes directly to
     * another.
     * @param thingReader the reader.
     * @param thingReceiver the receiver.
     */
    public ThingParser(ThingReader thingReader, ThingReceiver thingReceiver) {
        this(thingReader, thingReceiver, new ThingConverter<?>[0]);
    }

    /**
     * Construct a parser that reads from one source and writes to another
     * via a chain of converters.
     * @param thingReader the reader.
     * @param thingReceiver the receiver.
     * @param thingConverters the converters.
     */
    public ThingParser(ThingReader thingReader, ThingReceiver thingReceiver, ThingConverter<? extends Serializable>... thingConverters) {
        this.setReader(thingReader);
        this.setReceiver(thingReceiver);
        for (ThingConverter<? extends Serializable> thingConverter : thingConverters) {
            this.addConverter(thingConverter);
        }
    }
    
    /**
     * Calls the close() method on the reader, receiver, and all converters,
     * in order to tidy up after parsing and free resources.
     * @throws IOException if they could not be closed.
     */
    public void close() throws IOException {
        this.thingReader.close();
        for (ThingConverter<?> converter : this.thingConverters) {
            converter.close();
        }
        this.thingReceiver.close();
    }

    /**
     * Sets the reader to use as a source of things.
     * @param thingReader the reader.
     */
    public void setReader(ThingReader thingReader) {
        this.thingReader = thingReader;
    }

    /**
     * Sets the receiver to use to iterate over the resulting things. The 
     * parser will call 
     * {@link ThingParserMember#setNextThingReceiver(ThingReceiver)} on the 
     * reader or the last converter if it exists.
     * @param thingReceiver the receiver.
     * @exception IllegalStateException if the reader has not been set yet.
     */
    public void setReceiver(ThingReceiver thingReceiver) {
        if (this.thingReader == null) {
            throw new IllegalStateException("Cannot set receiver before reader.");
        }
        this.thingReceiver = thingReceiver;
        if (this.thingConverters.isEmpty()) {
            this.thingReader.setNextThingReceiver(thingReceiver);
        } else {
            this.thingConverters.get(this.thingConverters.size() - 1).setNextThingReceiver(thingReceiver);
        }
    }

    /**
     * Adds a converter to the end of the chain.
     * @param thingConverter the converter to add.
     * @exception IllegalStateException if the reader or receiver has not been
     * set yet.
     */
    public void addConverter(ThingConverter<? extends Serializable> thingConverter) {
        if (this.thingReader == null) {
            throw new IllegalStateException("Cannot set converter before reader.");
        }
        if (this.thingReceiver == null) {
            throw new IllegalStateException("Cannot set converter before receiver.");
        }
        this.thingConverters.add(thingConverter);
        thingConverter.setNextThingReceiver(this.thingReceiver);
        if (this.thingConverters.size() == 1) {
            this.thingReader.setNextThingReceiver(thingConverter);
        } else {
            this.thingConverters.get(this.thingConverters.size() - 2).setNextThingReceiver(thingConverter);
        }
    }

    /**
     * If {@link #hasNext()} returns false, or {@link #next()} returns 
     * {@code null}, check this first to see if there was an exception behind 
     * the problem.
     * @return {@code null} if there was no problem, otherwise returns the
     * last exception raised. This exception will remain in place until the 
     * next call to {@link #hasNext()} or {@link #next()}.
     */
    public Exception getLastException() {
        return this.lastException;
    }

    public boolean hasNext() {
        logger.fine("Checking to see if reader can continue reading.");
        this.lastException = null;
        try {
            return this.thingReader.canReadNextThing();
        } catch (final IOException ioe) {
            logger.log(Level.WARNING, "Reader threw an exception.", ioe);
            this.lastException = ioe;
        }
        logger.fine("Reader cannot continue reading.");
        return false;
    }

    @SuppressWarnings("unchecked")
	public T next() {
        logger.fine("Moving to next Thing from reader.");
        this.lastException = null;
        if (this.thingReceiver instanceof ThingBuilder<?>) {
            try {
                this.thingReader.readNextThing();
            } catch (final IOException ioe) {
                logger.log(Level.WARNING, "Reader threw an exception.", ioe);
                this.lastException = ioe;
            }
            logger.fine("Obtaining finished thing from reader.");
            ThingBuilder<T> builder = (ThingBuilder<T>)ThingBuilder.class.cast(this.thingReceiver);
            return builder.getFinishedThing();
        }
        throw new UnsupportedOperationException("Can't iterate if the receiver is not a builder.");
    }

    /**
     * Iterates through all the things available and pushes them through the
     * parser. 
     * @throws IOException directly from {@link ThingReader#canReadNextThing()}
     * and {@link ThingReader#readNextThing()}.
     */
    public void parseAll() throws IOException {
        while (this.thingReader.canReadNextThing()) {
            this.thingReader.readNextThing();
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Can't remove things from the parser.");
    }
}
