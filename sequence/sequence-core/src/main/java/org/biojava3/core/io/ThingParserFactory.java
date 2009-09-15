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
import java.util.Iterator;

/**
 * A convenience toolkit for constructing quick parsers. When passed a {@link ThingFormat}
 * it uses the classes returned by it in the way described by the {@link ThingFormat} interface
 * to construct simple parser chains. The returned parser can be iterated using {@link ThingParser#hasNext()}
 * and other methods from the {@link Iterator} interface if reading, or can be done all-at-once using
 * {@link ThingParser#parseAll()} if writing.
 * 
 * @author Richard Holland
 * @since 3.0
 */
public class ThingParserFactory {

    /**
     * Constructs a parser to read from a {@link File} in the given format.
     * @param <T> a format type, e.g. FASTA.
     * @param format the format to use, e.g. FASTA.format.
     * @param file the file to read.
     * @return the parser.
     */
    public static <T extends Serializable> ThingParser<T> getReadParser(ThingFormat<T> format, File file) {
        return _getReadParser(format, format.getFileReaderClass(), new Class[]{File.class}, new Object[]{file});
    }

    /**
     * Constructs a parser to read from an {@link InputStream} in the given format.
     * @param <T> a format type, e.g. FASTA.
     * @param format the format to use, e.g. FASTA.format.
     * @param is the input stream to read.
     * @return the parser.
     */
    public static <T extends Serializable> ThingParser<T> getReadParser(ThingFormat<T> format, InputStream is) {
        return _getReadParser(format, format.getFileReaderClass(), new Class[]{InputStream.class}, new Object[]{is});
    }

    /**
     * Constructs a parser to read from a {@link Reader} in the given format.
     * @param <T> a format type, e.g. FASTA.
     * @param format the format to use, e.g. FASTA.format.
     * @param rdr the reader to read.
     * @return the parser.
     */
    public static <T extends Serializable> ThingParser<T> getReadParser(ThingFormat<T> format, Reader rdr) {
        return _getReadParser(format, format.getFileReaderClass(), new Class[]{Reader.class}, new Object[]{rdr});
    }

    /**
     * Constructs a parser to write to a {@link File} in the given format.
     * @param <T> a format type, e.g. FASTA.
     * @param format the format to use, e.g. FASTA.format.
     * @param file the file to write to.
     * @param thing the object to write out.
     * @return the parser.
     */
    public static <T extends Serializable> ThingParser<T> getWriteParser(ThingFormat<T> format, File file, T thing) {
        return _getWriteParser(format, format.getFileWriterClass(), new Class[]{File.class}, new Object[]{file}, thing);
    }

    /**
     * Constructs a parser to write to an {@link OutputStream} in the given format.
     * @param <T> a format type, e.g. FASTA.
     * @param format the format to use, e.g. FASTA.format.
     * @param os the output stream to write to.
     * @param thing the object to write out.
     * @return the parser.
     */
    public static <T extends Serializable> ThingParser<T> getWriteParser(ThingFormat<T> format, OutputStream os, T thing) {
        return _getWriteParser(format, format.getFileWriterClass(), new Class[]{OutputStream.class}, new Object[]{os}, thing);
    }

    /**
     * Constructs a parser to write to a {@link Writer} in the given format.
     * @param <T> a format type, e.g. FASTA.
     * @param format the format to use, e.g. FASTA.format.
     * @param wtr the writer to write to.
     * @param thing the object to write out.
     * @return the parser.
     */
    public static <T extends Serializable> ThingParser<T> getWriteParser(ThingFormat<T> format, Writer wtr, T thing) {
        return _getWriteParser(format, format.getFileWriterClass(), new Class[]{Writer.class}, new Object[]{wtr}, thing);
    }

    /**
     * Build a generic parser which can read from any source.
     * @param <T> a format type, e.g. FASTA.
     * @param format the format to use, e.g. FASTA.format.
     * @param readerClass the class of reader to use.
     * @param paramTypes the parameter types of the reader class constructor.
     * @param params the params for the reader class constructor.
     * @return the parser.
     */
    private static <T extends Serializable> ThingParser<T> _getReadParser(ThingFormat<T> format, Class<? extends ThingReader> readerClass, Class<?>[] paramTypes, Object[] params) {
        Class<? extends ThingBuilder<T>> builderClass = format.getBuilderClass();
        try {
            ThingReader reader = (ThingReader) readerClass.getConstructor(paramTypes).newInstance(params);
            ThingBuilder<T> builder = builderClass.newInstance();
            return new ThingParser<T>(reader, builder);
        } catch (Exception e) {
            throw new IllegalArgumentException("Format does not comply with ThingParserFactory requirements", e);
        }
    }

    /**
     * Build a generic parser which can write to any source.
     * @param <T> a format type, e.g. FASTA.
     * @param format the format to use, e.g. FASTA.format.
     * @param writerClass the class of writer to use.
     * @param paramTypes the parameter types of the writer class constructor.
     * @param params the params for the writer class constructor.
     * @param thing the thing to write.
     * @return the parser.
     */
    private static <T extends Serializable> ThingParser<T> _getWriteParser(ThingFormat<T> format, Class<? extends ThingWriter> writerClass, Class<?>[] paramTypes, Object[] params, T thing) {
        Class<? extends ThingEmitter<T>> emitterClass = format.getEmitterClass();
        if (thing == null) {
            throw new IllegalArgumentException("Cannot write null objects");
        }
        try {
            ThingEmitter<T> emitter = emitterClass.getConstructor(new Class[]{thing.getClass()}).newInstance(new Object[]{thing});
            ThingWriter writer = writerClass.getConstructor(paramTypes).newInstance(params);
            return new ThingParser<T>(emitter, writer);
        } catch (Exception e) {
            throw new IllegalArgumentException("Format does not comply with ThingParserFactory requirements", e);
        }
    }
}
