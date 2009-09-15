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

import java.io.File;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Reader;
import java.io.Serializable;
import java.io.Writer;

/**
 * Each file format class should provide a static final constant variable in
 * the main format interface called 'format', which should be an instance
 * of a class implementing this interface. This should specify the default
 * classes to use when reading and writing format objects. When passed to
 * the {@link ThingParserFactory}, the responses will be used to construct
 * parser chains for simple read/write operations.
 * 
 * @author Richard Holland
 * @since 3.0
 */
public interface ThingFormat<T extends Serializable> {

    /**
     * Obtain the builder class. It should be able to operate having been
     * instantiated via a no-args constructor.
     * @return the builder class.
     */
    public Class<? extends ThingBuilder<T>> getBuilderClass();

    /**
     * Obtain the emitter class. It should be able to operate having been
     * instantiated via a constructor taking just one arg, which is an 
     * object of the file format this belongs to (T).
     * @return the emitter class.
     */
    public Class<? extends ThingEmitter<T>> getEmitterClass();

    /**
     * Obtain the file reader class. It should be able to operate having been
     * instantiated via a constructor taking just one arg, which is a {@link File}.
     * @return the reader class.
     */
    public Class<? extends ThingReader> getFileReaderClass();

    /**
     * Obtain the stream reader class. It should be able to operate having been
     * instantiated via a constructor taking just one arg, which is an {@link InputStream}.
     * @return the reader class.
     */
    public Class<? extends ThingReader> getStreamReaderClass();

    /**
     * Obtain the reader reader class. It should be able to operate having been
     * instantiated via a constructor taking just one arg, which is a {@link Reader}.
     * @return the reader class.
     */
    public Class<? extends ThingReader> getReaderReaderClass();

    /**
     * Obtain the file writer class. It should be able to operate having been
     * instantiated via a constructor taking just one arg, which is a {@link File}.
     * @return the writer class.
     */
    public Class<? extends ThingWriter> getFileWriterClass();

    /**
     * Obtain the stream writer class. It should be able to operate having been
     * instantiated via a constructor taking just one arg, which is a {@link OutputStream}.
     * @return the writer class.
     */
    public Class<? extends ThingWriter> getStreamWriterClass();

    /**
     * Obtain the writer writer class. It should be able to operate having been
     * instantiated via a constructor taking just one arg, which is a {@link Writer}.
     * @return the writer class.
     */
    public Class<? extends ThingWriter> getWriterWriterClass();
}
