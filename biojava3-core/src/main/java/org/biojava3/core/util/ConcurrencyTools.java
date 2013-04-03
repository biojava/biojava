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
 * Created on May 26, 2010
 * Author: Mark Chapman
 */

package org.biojava3.core.util;

import java.util.concurrent.Callable;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/**
 * Static utility to easily share a thread pool for concurrent/parallel/lazy execution.  To exit cleanly,
 * {@link #shutdown()} or {@link #shutdownAndAwaitTermination()} must be called after all tasks have been submitted.
 *
 * @author Mark Chapman
 */
public class ConcurrencyTools {

    private static ThreadPoolExecutor pool;
    // private static int tasks;

    // TODO additional logging and listening services

    // prevents instantiation
    private ConcurrencyTools() { }

    /**
     * Returns current shared thread pool.  Starts up a new pool, if necessary.
     *
     * @return shared thread pool
     */
    public static ThreadPoolExecutor getThreadPool() {
        if (pool == null || pool.isShutdown()) {
            setThreadPoolDefault();
        }
        return pool;
    }

    /**
     * Sets thread pool to reserve a given number of processor cores for foreground or other use.
     *
     * @param cpus number of processor cores to reserve
     */
    public static void setThreadPoolCPUsAvailable(int cpus) {
        setThreadPoolSize(Math.max(1, Runtime.getRuntime().availableProcessors() - cpus));
    }

    /**
     * Sets thread pool to a given fraction of the available processors.
     *
     * @param fraction portion of available processors to use in thread pool
     */
    public static void setThreadPoolCPUsFraction(float fraction) {
        setThreadPoolSize(Math.max(1, Math.round(fraction * Runtime.getRuntime().availableProcessors())));
    }

    /**
     * Sets thread pool to default of 1 background thread for each processor core.
     */
    public static void setThreadPoolDefault() {
        setThreadPoolCPUsAvailable(0);
    }

    /**
     * Sets thread pool to a single background thread.
     */
    public static void setThreadPoolSingle() {
        setThreadPoolSize(1);
    }

    /**
     * Sets thread pool to given size.
     *
     * @param threads number of threads in pool
     */
    public static void setThreadPoolSize(int threads) {
    	setThreadPool(   new ThreadPoolExecutor(threads, threads,
                                      0L, TimeUnit.MILLISECONDS,
                                      new LinkedBlockingQueue<Runnable>()));
    	
        
    }

    /**
     * Sets thread pool to any given {@link ThreadPoolExecutor} to allow use of an alternative execution style.
     *
     * @param pool thread pool to share
     */
    public static void setThreadPool(ThreadPoolExecutor pool) {
        if (ConcurrencyTools.pool != pool) {
            shutdown();
            ConcurrencyTools.pool = pool;
        }
    }

    /**
     * Disables new tasks from being submitted and closes the thread pool cleanly.
     */
    public static void shutdown() {
        if (pool != null) {
            pool.shutdown();
        }
    }

    /**
     * Closes the thread pool.  Waits 1 minute for a clean exit; if necessary, waits another minute for cancellation.
     */
    public static void shutdownAndAwaitTermination() {
        shutdown();
        if (pool != null) {
            try {
                // wait a while for existing tasks to terminate
                if (!pool.awaitTermination(60L, TimeUnit.SECONDS)) {
                    pool.shutdownNow(); // cancel currently executing tasks
                    // wait a while for tasks to respond to being canceled
                    if (!pool.awaitTermination(60L, TimeUnit.SECONDS)) {
                        System.err.println("BioJava ConcurrencyTools thread pool did not terminate");
                    }
                }
            } catch (InterruptedException ie) {
                pool.shutdownNow(); // (re-)cancel if current thread also interrupted
                Thread.currentThread().interrupt(); // preserve interrupt status
            }
        }
    }

    /**
     * Queues up a task and adds a log entry.
     *
     * @param <T> type returned from the submitted task
     * @param task submitted task
     * @param message logged message
     * @return future on which the desired value is retrieved by calling get()
     */
    public static<T> Future<T> submit(Callable<T> task, String message) {
        // TODO log("Task " + (++tasks) + " submitted to shared thread pool. " + message);
        return getThreadPool().submit(task);
    }

    /**
     * Queues up a task and adds a default log entry.
     *
     * @param <T> type returned from the submitted task
     * @param task submitted task
     * @return future on which the desired value is retrieved by calling get()
     */
    public static<T> Future<T> submit(Callable<T> task) {
        return submit(task, "");
    }

}
