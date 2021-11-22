TCO 21 Final problem: CoinCollector
Final Rank: 1

*tl;dr* I added my commented code here: https://github.com/CatalinTiseanu/topcoder-marathon/blob/main/tco21-final/CoinCollector.cpp
It's fairly well commented and starts at line ~550.
The main action happens in Rec3() and in Solve().

# Problem & thought process

After 30+ rounds at 5 prior TCO finals this is definitely the most special 'Post your approach post' I wrote so far :)

Fun fact: When I first see a problem statement I check to see if it's deterministic or has randomness and if it has randomness I semi-panic. I totally bombed on them in 2018 missing the qualification to the final. I guess this time around the randomness was more 'predictable' and didn't impact the LB quite that much.

I really liked that this problem was really easy to understand from the get go; I barely read the visualizer code just to make sure I didn't misunderstand anything. Getting something going was really easy this time around and this problem could have worked as well for a 10 hour final.

I first thought about the possibility of beam search, but it felt that since you have so many replacement possibilities and randomness it might not be so effective.
I also felt that 'exact' computation for a few turns in the future is more important than imperfect simulations for a few more turns (even though I didn't try this).
Finally, when I saw the possible combinations of dice I thought about some DP around (dice in hand, cell) however I didn't see how it would help with the coins. I should have revisited this towards the end (I talk more about this in Improvements). 

# What I did in 24 hours (felt like so much time and I never wanted a round to finish earlier that bad)

I started with the lookahead idea and kept at it. The final started at 11 am local time here and I managed to grab some sleep round midnight for 4 hours.
The fact that I was able to fall asleep with everything that was going on was really surprising and a secondary achievement for this final :)
It was really hard to keep motivated, but knowing you guys I was paranoid until the end. Psyho really made it tough at the end since due to variance 2 of my last solutions were worse than him.
My final solution only used 9 seconds, even though there was a guaranteed anti-TLE strategy in just finishing the simulation.
However since I couldn't test on a c3 and the artifact logs were invalid (not founds) for my last submissions I decided that I'd rather lose by a small margin than due to TLE. using an addition 0.7 seconds would have helped for sure.


# Solution

The crux of my solution was a depth (lookhead) search.
Intutively you can think of this as two player game where the adversary chooses the value of the dice randomly.
The search then prioritized solely getting the maximum possible score over the search window (which was 4 at most).
No long-term calculation was made for leafs.
As in alpha-beta style two-player games, pruning early was key.
Finally, all the tuning I did only looked at the number of dices (no tuning based on % coins, % spikes or N).


The main optimizations were around:
1. What replacement dice to consider after every recursive call
2. How to estimate when it's worth continuing after hitting a spike (when the expected score at the end of the search would be less than the current score prior to making any move)
3. What small bonuses / maluses to add in order to break tiebreaks in the leaf evaluation
4. How to select lookahead depth and adjust this depending on remaining time
5. How to remove parts of the search tree efficiently (prune)


1. I always use the STAY dice when I get it unless it's towards the end of the game (in the last few turns it might be better to try to get that last coin instead of wasting a turn with STAY although I didn't test this).
I always assume I get RANDOM as replacement except when D <= 3 when for the very first replacement I try out all 12 possibilities before continuing the search. This gave me a 6-7% boost for 2 dices and and based on [handle]wleite[/handle] stats (thank you btw!) it might have been one of the decisive optimizations.

2. I compute how much score per turn I had so far (current_score / # of turns).
However this alone is not a good estimate for the future since collecting coins becomes harder and harder as the turns go on.
Therefore, I need to adjust this score to account for that, for which I use fudge_factor which is
```
fudgefactor = 0.3 + 0.1 * D;
        if (D == 2)
            fudgefactor = 0.3;
```    
        
I also compute max remaining score which depends on the number of remaining coins as well as number of turns.

I need that the expected_score + fudge_factor * remaining_max_score >= current_score in order to continue the search.
This is basically what controls whether I continue after hitting a spike.

3. If D >= 4 I give a bonus if using the RANDOM dice (reducing randomness) as well as a malus for hitting a spike

4. I start with lookahead = 4, and then compute expected remaining time (based # of turns = A -  current_turn) and if it doesn't look I switch to lookahead = 3. Once in panic mode lookahead = 2. I stressed out a lot around TLE since the instances we were allowed to use were different than the judge one (which was a c3). Note for future competitions: make sure the instance the competitors can use is the same as the judge one. 
A trivial fix (which believe or not I thought about only 10 mins before finish LOL) was to just exist the simulation once time was close 

5. One of the key additions was around handling spikes. I allow moves into spikes however I first compute the maximum possible score. To give an example, let's say I want to go RIGHT with a RANDOM dice and 3 our of the 6 possible landings have a spike. Let's say I have 2 more steps in my search and I have enough coins around. I compute the maximum expected score as 
(3/6) * (current_score + 1 * 100) // I hit a spike and right after I still collect a coin
+ (3/6) * (current_score + 2 * 100) // I don't hit a spike and get max score.

Looking at this now after writing it I could obviously take into account whether the non-spike cell I land actually has a coin or not instead of just assuming it has but oh well :)

The key thing is that this estimation is really quick to compute and I use it to prune a lot of the search space.

# Improvements and closing thoughts

A big weakness of my approach is that there is no sense of a long-term value for a cell (let's say a leaf in my search tree).
I want to revisit the DP approach to compute at the beginning, taking just spikes into account, a safety measure of the form
'expected number of steps from this state until hitting a spike' and add this tl 

Thank you to the organizers for another year of great problems; I really like the return (valid qualification problems as well) to more visual style problems where the visualizer is important and helps you get some intuition around what is happening :)

To many more TCO Marathon onsite finals!
