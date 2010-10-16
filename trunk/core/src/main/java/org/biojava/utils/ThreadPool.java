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

package org.biojava.utils;

/**
 * <p><code>ThreadPool</code> specifies basic thread-pooling
 * operations such that third-party implementations may be used
 * without requiring changes to BioJava.</p>
 *
 * @author Keith James
 * @since 1.3
 */
public interface ThreadPool
{
    /**
     * <code>addRequest</code> requests that a <code>Runnable</code>
     * be scheduled to be run by one of the threads in the pool.
     *
     * @param task a <code>Runnable</code>.
     */
    public void addRequest(Runnable task);

    /**
     * <code>startThreads</code> starts all the threads running and
     * opens the pool to requests.
     */
    public void startThreads();

    /**
     * <code>stopThreads</code> causes all running threads to stop
     * when their current task is complete. It also closes the pool to
     * new requests. Requests still queued are not done and the queue
     * is emptied.
     */
    public void stopThreads();

    /**
     * <code>waitForThreads</code> temporarily closes the pool to new
     * requests until such time as the current request queue has been
     * emptied and all running tasks completed.
     */
    public void waitForThreads();
}
