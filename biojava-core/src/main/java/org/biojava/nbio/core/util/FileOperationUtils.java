package org.biojava.nbio.core.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;

/**
 * @author Sreejith Nair
 */
public class FileOperationUtils {
    /**
     * @param src
     * @param dst
     * @throws IOException
     */
    @SuppressWarnings("resource")
    public static void copy(File src, File dst) throws IOException {

        // Took following recipe from
        // http://stackoverflow.com/questions/106770/standard-concise-way-to-copy-a-file-in-java
        // The nio package seems to be the most efficient way to copy a file
        FileChannel source = null;
        FileChannel destination = null;

        try {
            // we need the supress warnings here (the warning that the stream is not closed is harmless)
            // see http://stackoverflow.com/questions/12970407/does-filechannel-close-close-the-underlying-stream
            source = new FileInputStream(src).getChannel();
            destination = new FileOutputStream(dst).getChannel();
            destination.transferFrom(source, 0, source.size());
        } finally {
            if (source != null) {
                source.close();
            }
            if (destination != null) {
                destination.close();
            }
        }
    }

    /**
     * Gets the file extension of a file, excluding '.'.
     * If the file name has no extension the file name is returned.
     * @param f a File
     * @return The extension
     */
    public static String getFileExtension(File f) {
        String fileName = f.getName();
        String ext = "";
        int mid = fileName.lastIndexOf(".");
        ext = fileName.substring(mid + 1, fileName.length());
        return ext;
    }

    /**
     * Gets the file name up to and excluding the first
     * '.' character. If there is no extension, the full filename
     * is returned.
     * @param f A file
     * @return A possibly empty but non-null String.
     */
    public static String getFilePrefix(File f) {
        String fileName = f.getName();
        int mid = fileName.indexOf(".");
        if (mid < 0) {
            return fileName;
        }
        return fileName.substring(0, mid);
    }

    /**
     * Converts path to Unix convention and adds a terminating slash if it was
     * omitted.
     *
     * @param path original platform dependent path
     * @return path in Unix convention
     * @author Peter Rose
     * @since 3.2
     */
    public static String toUnixPath(String path) {
        String uPath = path;
        if (uPath.contains("\\")) {
            uPath = uPath.replaceAll("\\\\", "/");
        }
        // this should be removed, it's need since "\" is added AtomCache code
        if (uPath.endsWith("//")) {
            uPath = uPath.substring(0, uPath.length() - 1);
        }
        if (!uPath.endsWith("/")) {
            uPath = uPath + "/";
        }
        return uPath;
    }

    /**
     * Expands ~ in paths to the user's home directory.
     *
     * <p>
     * This does not work for some special cases for paths: Other users' homes
     * (~user/...), and Tilde expansion within the path (/.../~/...). In these cases
     *  the original argument is returned.
     *
     * @param file A filepath starting with a tilde
     * @return An absolute path
     */
    public static String expandUserHome(String file) {
        // replace any / with the proper separator (/ or \ for Linux and Windows respectively).
        file = file.replaceAll("/", "\\"+File.separator); //The "\\" is to escape the separator if needed.
        if (file.startsWith("~") && (file.length() == 1 || File.separator.equals(file.substring(1, 2)))) {
            file = System.getProperty("user.home") + file.substring(1);
        }
        return file;
    }
    /**
     * Recursively delete a folder & contents
     *
     * @param dir directory to delete
     */
    public static void deleteDirectory(Path dir) throws IOException {
        if(dir == null || !Files.exists(dir))
            return;
        Files.walkFileTree(dir, new SimpleFileVisitor<Path>() {
            @Override
            public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                Files.delete(file);
                return FileVisitResult.CONTINUE;
            }

            @Override
            public FileVisitResult postVisitDirectory(Path dir, IOException e) throws IOException {
                if (e != null) {
                    throw e;
                }
                Files.delete(dir);
                return FileVisitResult.CONTINUE;
            }
        });
    }
    /**
     * Recursively delete a folder & contents
     *
     * @param dir directory to delete
     */
    public static void deleteDirectory(String dir) throws IOException {
        deleteDirectory(Paths.get(dir));
    }
}
