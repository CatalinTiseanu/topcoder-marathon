TCO 21 Final problem: CoinCollector Final Rank: 1
**tl;dr** I added my commented code here: https://github.com/CatalinTiseanu/topcoder-marathon/blob/main/tco21-final/CoinCollector.cpp It's reasonably well commented and starts at line ~550. The main action happens in Rec3() and in Solve().

After 30+ rounds at five prior TCO finals, this is the most memorable 'Post your approach post' I have written so far :)

# Problem & thought process

**Fun fact:** When I first see a problem statement, I check to see if it's deterministic or has randomness - if it does, I get scared. I bombed on the randomness rounds in 2018, missing the qualification to the final. This time around, randomness was more 'predictable' and didn't impact the LB quite that much.
I liked that this problem was straightforward to understand from the get-go; I barely read the visualizer code to ensure I didn't misunderstand anything. Getting something going was easy this time around, and this problem could have worked as well for a 10-hour final.

I first thought about the possibility of beam search, but I felt that it might not be so effective since you have so many replacement possibilities and randomness. I also felt that 'exact' computation for a few turns in the future is more critical than imperfect simulations for a few more turns (even though I didn't try this). 

Finally, when I saw the possible dice combinations, I thought about some DP around (dice in hand, cell); however, I didn't see how it would help with the coins. I should have revisited this towards the end (I talk more about this in Improvements).

# What I did in 24 hours (felt like so much time, and I never wanted a contest to finish earlier that bad)

I started with the lookahead idea and kept at it. The final started at 11 am local time here, and I managed to grab some sleep around midnight for 4 hours. The fact that I was able to fall asleep with everything that was going on was surprising and a secondary achievement for this final :) It was tough to keep motivated, but knowing you guys, I was paranoid until the end. Psyho made it challenging at the end since, due to variance, two of my last solutions were worse than him. My final solution only used 9 seconds, even though there was a guaranteed anti-TLE strategy in just finishing the simulation. However, since I couldn't test on a c3 and the artifact logs were invalid ('not found') for my last submissions, I decided that I'd rather lose by a small margin than due to TLE. Using an additional 0.7 seconds would have helped for sure.

# Solution

The crux of my solution was a depth (lookahead) search. Intuitively, you can think of this as a two-player game where the adversary randomly chooses the dice's value. The search then prioritized solely getting the maximum possible score over the search window (which was four at most). I made no long-term calculation for leaves. As in 'alpha-beta'-style two-player games, pruning early was important. Finally, all the tuning I did only looked at the number of dices (no tuning based on % coins, % spikes, or N).

The main optimizations were around:

1) What replacement dice to consider after every recursive call
2 ) How to estimate when it's worth continuing after hitting a spike (when the expected score at the end of the search would be less than the current score before making any move)
3) What small bonuses/maluses to add to break tiebreaks in the leaf evaluation
4) How to select lookahead depth and adjust this depending on the remaining time
5) How to remove parts of the search tree efficiently (prune)

## My answers were:

1) I always use the STAY dice when I get it unless it's towards the end of the game (in the last few turns, it might be better to try to get that final coin instead of wasting a turn with STAY although I didn't test this). I always assume I get RANDOM as a replacement except when D <= 3 when for the very first replacement, I try out all 12 possibilities before continuing the search. This extra logic gave me a 6-7% boost for two dices and based on [handle]wleite[/handle] stats (thank you, btw!), it might have been one of the decisive optimizations.

2) I compute how much score per turn I had so far (current_score / # of turns). However, this alone is not a decent estimate for the future since collecting coins becomes harder and more complex as the turns go on. Therefore, I need to adjust this score to account for that, for which I use fudge_factor, which in code is
```
	fudgefactor = 0.3 + 0.1 * D;
        if (D == 2)
            fudgefactor = 0.3;
 ```
I also compute the max remaining score, which depends on the number of remaining coins and the number of turns.
I need that the expected_score + fudge_factor * remaining_max_score >= current_score in order to continue the search. This formula is what controls whether I continue after hitting a spike.

3) If D >= 4, I give a bonus if using the RANDOM dice (reducing randomness) as well as a malus for hitting a spike

4) I start with lookahead = 4 and then compute the expected remaining time (based # of turns = A - current_turn), and if it doesn't look like there will be enough time, I switch to lookahead = 3. Once in panic mode lookahead becomes 2. I stressed out a lot around TLE since the instances we were allowed to use were different from the judge one (a c3). Note for future competitions: make sure the competitors' instance is the same as the judge one. A trivial fix (which, believe or not, I thought about only 10 mins before finishing - LOL) was to exit the simulation once the time was close.

5) One of the key additions was around handling spikes. I allow moves into spikes; however, I first compute the maximum possible score. To give an example, let's say I want to go RIGHT with a RANDOM dice and 3 our of the 6 possible landings have a spike. Let's say I have 2 more steps in my search, and I have enough coins around. I compute the maximum expected to score as the sum of
* (3/6) * (current_score + 1 * 100) // I hit a spike and right after I still collect a coin
* (3/6) * (current_score + 2 * 100) // I don't hit a spike and get max score.
Looking at this now after writing it, I could take into account whether the non-spike cell I land actually has a coin instead of just assuming it has, but oh well :)
The critical thing is that this estimation is rapid to compute, and I use it to prune the search space a lot.

# Improvements and closing thoughts

A significant weakness of my approach is that there is no sense of a long-term value for a cell (let's say a leaf in my search tree). I want to revisit the DP approach to compute at the beginning, taking just spikes into account, a safety measure of the form 'expected number of steps from this state until hitting a spike' and add this to the evaluation score.

**Thank you to the organizers for another year of interesting problems; I really enjoy the return (valid qualification problems as well) to more visual style problems where the visualizer is essential and helps you get some intuition around what is happening :) 

To many more TCO Marathon onsite finals!**
