import argparse
import sys
import math
from gzip import GzipFile



# Prepare for cmd line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-f", "--file", help="Input simulation trace file", default=None)
parser.add_argument("-z", "--zip", help="Read input trace from a zipped file", default=False, action='store_const', const=True)
args = parser.parse_args()

# Are we reading from a file specified by the user, or from the standard input?
if args.file != None:

    try:

        f = open(args.file, "r")
    except IOError:

        print ("Could not open file '{0}'".format(args.file))
else:

    f = sys.stdin

if args.zip == True:
    f = GzipFile(fileobj=f, mode="r")

# Initialize some data
genPackets = {}
mediumAccessPackets = {}
deliveredPackets = {}
authorizedPackets = {}
abortedPackets = {}
senderSuccessful = {}
senderDrops = {}
senderAttempt = {}
receiverAckSuccess = {}
receiverSuccess = {}
ackAttempt = {}
totalAttempts = 0
totalAckAttempts = 0
receivedWithCollision = 0
droppedWithCollision = 0
receivedWithoutCollision = 0
# Main loop: iterate through the lines of the input.
data = f.readlines()
for line in data:

    # Each entry in the log is an event with a well-defined format and fields
    # separated by spaces. Let's start by splitting those fields.
    items = line.split()

    # First step of processing this event: find out its type.
    type = items[0]
    if type == "+":

        # New packet generated at one of the stations. Generate when that
        # happened.
        genPackets[items[2] + items[3]] = float(items[1])
    elif type == "G":

        # The process of transmitting an application layer packet begins now for
        # the link layer. Store this time for delay-computing purposes.
        authorizedPackets[items[2] + items[3]] = float(items[1])
    elif type == "A":

        # Transmission of the application level packet has been aborted. Store
        # this information for packet loss related computation.
        abortedPackets[items[2] + items[3]] = float(items[1])
    elif type == "S":

        # Sender understands that this application layer packet has been successfully
        # delivered to the AP.
        senderSuccessful[items[2] + items[3]] = float(items[1])
    elif type == "D":

        # Sender thinks this packet has not been delivered to destination and is
        # giving up.
        senderDrops[items[2] + items[3]] = float(items[1])
    elif type == "MDs":

        # Start of a DIFS wait for this packet. Let's store the time this
        # event happened for the first time for each packet.
        if not items[2] + items[3] in mediumAccessPackets:
            mediumAccessPackets[items[2] + items[3]] = float(items[1])
            #print("oiiiiiiiiiii")

    elif type == "To" and len(items) == 4:

        # Link-layer transmission attempt for this packet
        if not items[2] + items[3] in senderAttempt:
            senderAttempt[items[2] + items[3]] = []

        senderAttempt[items[2] + items[3]].append(float(items[1]))
        totalAttempts = totalAttempts + 1
    elif type == "To":

        # Link-layer transmission attempt for an ack
        if not items[3] + items[4] in ackAttempt:
            ackAttempt[items[3] + items[4]] = []

        ackAttempt[items[3] + items[4]].append(float(items[1]))
        totalAckAttempts = totalAckAttempts + 1
    elif type == "r" and items[4] != "[ack]":

        # Link-layer data transmission attempt was successful.
        # Store only the first instance for each packet.
        if not items[3] + items[4] in receiverSuccess:
            receiverSuccess[items[3] + items[4]] = float(items[1])

        if int(items[5]) > 1:
            receivedWithCollision = receivedWithCollision + 1
        else:
            receivedWithoutCollision = receivedWithoutCollision + 1
    elif type == "r":

        # Link-layer ack transmission attempt was successful.
        receiverAckSuccess[items[2] + items[3]] = float(items[1])
    elif type == "d" and items[4] != "[ack]":

        # Link-layer data transmission attempt was unsuccessful.
        if int(items[5]) > 1:
            droppedWithCollision = droppedWithCollision + 1


# After processing all events on the log, compute desired metrics.
print ("#####################")
print ("Losses and deliveries")
print ("#####################")
print ("### Application layer statistics:\n")
nGenPackets = len(genPackets)
print ("\tTotal number of application layer packets: {0}".format(nGenPackets))

nAborted = len(abortedPackets)
print ("\tNumber of packets that were eventually aborted: {0}".format(nAborted))

nMediumAcess = len(mediumAccessPackets)
print ("\tNumber of packets that eventually reached the medium access phase (but could have been aborted later): {0}".format(nMediumAcess))

nAppReceived = len(receiverSuccess)
print ("\tNumber of packets that were actually delivered at the receiver's application layer: {0}".format(nAppReceived))

print ("\tPacket delivery rate wrt total number of generated packets: {0}".format(float(nAppReceived) / nGenPackets))
print ("\tPacket delivery rate wrt those that actually reached medium access: {0}".format(float(nAppReceived) / (nMediumAcess)))

print ("\n### Link-layer statistics:\n")
print ("\tData frame transmission attempts: {0}".format(totalAttempts))
print ("\tData frame transmission successes (only forward direction): {0}".format(totalAckAttempts))
print ("\tData frame ack receptions (complete success): {0}".format(len(receiverAckSuccess)))
print ("\tFrame delivery rate (forward): {0}".format(totalAckAttempts / float(totalAttempts)))
print(totalAckAttempts)
print ("\tFrame delivery rate (backward): {0}".format(len(receiverAckSuccess) / float(totalAckAttempts)))
print ("\tFrame delivery rate (bidirectional): {0}".format(len(receiverAckSuccess) / float(totalAttempts)))

print ("\n#####################")
print ("Delay")
print ("#####################")
print ("### Application layer statistics:\n")

delaySum = 0
delaySquaredSum = 0
for i in receiverSuccess:
    sample = receiverSuccess[i] - genPackets[i]
    delaySum = delaySum + sample
    delaySquaredSum = delaySquaredSum + math.pow(sample, 2)

avgDelay = delaySum / len(receiverSuccess)
avgDelaySquared = delaySquaredSum / len(receiverSuccess)
varDelay = avgDelaySquared - avgDelay
stdevDelay = math.sqrt(varDelay)

print ("Average: {0} us".format(avgDelay))
print ("Standard deviation: {0} us".format(stdevDelay))

print ("\n#####################")
print ("Collisions")
print ("#####################")
print ("Number of data frames received withOUT collision: {0}".format(receivedWithoutCollision))
print ("Number of data frames received even with collision: {0}".format(receivedWithCollision))
print ("Number of data frames dropped with collision: {0}".format(droppedWithCollision))
print ("Total fraction of collided data frames: {0}".format((float(droppedWithCollision + receivedWithCollision)) / totalAttempts))
