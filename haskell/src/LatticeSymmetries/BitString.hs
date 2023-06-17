-- |
-- Module      : LatticeSymmetries.BitString
-- Description : Defines the 'BitString' data type
-- Copyright   : (c) Tom Westerhout, 2022
-- Stability   : experimental
module LatticeSymmetries.BitString
  ( BitString (..),
    readBitString,
    writeBitString,
    writeManyBitStrings,
    wordsFromBitString,
  )
where

import Control.Monad.Primitive
import Control.Monad.ST
import Data.Bits
import qualified Data.Primitive.Ptr as P
import Data.Vector.Generic (Vector)
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import Foreign.Ptr

-- | A newtype over 'Integer' indicating that we don't treat it as a number, but rather as a
-- collection of bits.
--
-- The order of bits is determined by the 'Bits' instance of 'Integer'.
newtype BitString = BitString {unBitString :: Integer}
  deriving stock (Show, Eq, Ord)
  deriving newtype (Bits)

writeBitString :: PrimMonad m => Int -> Ptr Word64 -> BitString -> m ()
writeBitString numberWords p (BitString x₀) = do
  let go !i !x
        | i < numberWords = do
          P.writeOffPtr p i (fromIntegral x :: Word64)
          go (i + 1) (x `shiftR` 64)
        | otherwise = pure ()
  go 0 x₀

writeManyBitStrings :: PrimMonad m => Int -> Ptr Word64 -> [BitString] -> m ()
writeManyBitStrings numberWords p xs = go p xs
  where
    go !ptr (y : ys) = writeBitString numberWords ptr y >> go (P.advancePtr ptr numberWords) ys
    go _ [] = pure ()

readBitString :: PrimMonad m => Int -> Ptr Word64 -> m BitString
readBitString numberWords p = do
  let go !i !x
        | i < numberWords = do
          y <- BitString <$> (fromIntegral :: Word64 -> Integer) <$> P.readOffPtr p i
          go (i + 1) ((y `shiftL` (64 * i)) .|. x)
        | otherwise = pure x
  go 0 zeroBits

wordsFromBitString :: Vector v Word64 => Int -> BitString -> v Word64
wordsFromBitString numberWords (BitString x₀) = runST $ do
  mvector <- GM.unsafeNew numberWords
  let go !i !x
        | i < numberWords = do
          GM.unsafeWrite mvector i (fromIntegral x :: Word64)
          go (i + 1) (x `shiftR` 64)
        | otherwise = pure ()
  go 0 x₀
  G.unsafeFreeze mvector
